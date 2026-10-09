"""Make the ``src`` layout importable without an editable install, and keep matplotlib headless."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

matplotlib.use("Agg", force=True)
