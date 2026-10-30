"""Hugging Face Space entry point (Gradio SDK).

The Space runs `python app.py` and serves whatever listens on port 7860. That is the
viraprime web app itself (FastAPI); MAFFT comes from packages.txt and the Python libraries
from requirements.txt. The `viraprime` package folder is uploaded next to this file.
"""

import os

import uvicorn

os.environ.setdefault("VIRAPRIME_JOBS_DIR", "/tmp/viraprime-jobs")

from viraprime.web.app import app  # noqa: E402  (after the environment is set)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("PORT", "7860")))
