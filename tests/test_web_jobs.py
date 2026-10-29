from __future__ import annotations

import json
import os
import time
from pathlib import Path

from viraprime.web.jobs import JobManager
from viraprime.web.params import Limits


def fake_job(job_dir: Path, argv: list[str]) -> None:
    """Sleeps argv[0] seconds, then writes a result like run_job does."""
    time.sleep(float(argv[0]))
    (Path(job_dir) / "result.json").write_text(json.dumps(
        {"exit_code": int(argv[1]) if len(argv) > 1 else 0, "mode": "tile", "rows": [],
         "error": "bad input" if len(argv) > 1 and argv[1] == "1" else None}))


def _wait(mgr: JobManager, job_id: str, states: set[str], timeout: float = 30) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        mgr.tick()
        st = mgr.status(job_id)
        if st["status"] in states:
            return st
        time.sleep(0.05)
    raise AssertionError(f"job {job_id} stuck in {mgr.status(job_id)}")


def test_lifecycle(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=1), target=fake_job)
    job = mgr.submit({"input.fasta": b">a\nACGT\n"}, ["0.5"])
    assert len(job) == 32 and (tmp_path / job / "input.fasta").read_bytes() == b">a\nACGT\n"
    assert mgr.status(job)["status"] == "queued"
    assert _wait(mgr, job, {"running", "done"})["status"] in {"running", "done"}
    st = _wait(mgr, job, {"done"})
    assert st["result"]["mode"] == "tile" and st["error"] is None


def test_exit_code_1_is_failed(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=1), target=fake_job)
    job = mgr.submit({}, ["0", "1"])
    st = _wait(mgr, job, {"done", "failed"})
    assert st["status"] == "failed" and st["error"] == "bad input"


def test_two_workers_and_queue_position(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=2), target=fake_job)
    jobs = [mgr.submit({}, ["3"]) for _ in range(3)]
    mgr.tick()
    assert [mgr.status(j)["status"] for j in jobs] == ["running", "running", "queued"]
    assert mgr.status(jobs[2])["queue_position"] == 1
    assert _wait(mgr, jobs[2], {"done"}, timeout=60)["status"] == "done"


def test_timeout(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=1, timeout_s=1), target=fake_job)
    job = mgr.submit({}, ["30"])
    st = _wait(mgr, job, {"failed", "done"})
    assert st["status"] == "failed" and "took longer than" in st["error"]


def test_progress_is_reported(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=1), target=fake_job)
    job = mgr.submit({}, ["3"])
    _wait(mgr, job, {"running"})
    (tmp_path / job / "progress.json").write_text('{"stage": "scoring", "progress": 0.4}')
    st = mgr.status(job)
    assert (st["stage"], st["progress"]) == ("scoring", 0.4)


def test_cleanup(tmp_path: Path) -> None:
    mgr = JobManager(tmp_path, Limits(workers=1, retention_s=100), target=fake_job)
    job = mgr.submit({}, ["0"])
    _wait(mgr, job, {"done"})
    stray = tmp_path / ("f" * 32)
    stray.mkdir()
    old = time.time() - 1000
    os.utime(stray, (old, old))
    assert mgr.cleanup() == 1
    assert not stray.exists() and (tmp_path / job).exists()
    assert mgr.cleanup(now=time.time() + 1000) == 1
    assert not (tmp_path / job).exists()
    try:
        mgr.status(job)
    except KeyError:
        pass
    else:
        raise AssertionError("expired job should be unknown")
