"""Reading measurements into the canonical in-memory intensity stack.

The file on disk is never modified, normalized, clipped or background-corrected.
Frames are returned as ``(N, H, W)`` in the raster LED order of the file.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy.io


@dataclass(frozen=True)
class IntensityData:
    values: np.ndarray          # (N, H, W) intensity, raster LED order
    source_path: Path | None


def _load_array(path: Path, key: str) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix == ".mat":
        content = scipy.io.loadmat(path)
        if key not in content:
            raise KeyError(f"MAT key {key!r} not found in {path}")
        return np.asarray(content[key])
    if suffix == ".npy":
        return np.asarray(np.load(path, allow_pickle=False))
    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as content:
            if key not in content:
                raise KeyError(f"NPZ key {key!r} not found in {path}")
            return np.asarray(content[key])
    raise ValueError(f"Unsupported data format {suffix!r}: {path}")


def load_intensity(config: dict[str, Any]) -> IntensityData:
    """Read ``config.data.file`` and reorder its axes to ``(N, H, W)``."""
    data = config["data"]
    path = Path(config["dataset_dir"]) / data["file"]
    if not path.is_file():
        raise FileNotFoundError(f"Measurement file not found: {path}")
    values = _load_array(path, str(data["key"]))

    layout = str(data["layout"]).upper()
    if values.ndim != 3:
        raise ValueError(f"Expected a 3-D stack in {path}, got shape {values.shape}")
    values = np.transpose(values, [layout.index(axis) for axis in "IHW"])
    if values.shape[1] != values.shape[2]:
        raise ValueError(f"Frames must be square, got {values.shape[1:]} in {path}")
    return IntensityData(values=np.ascontiguousarray(values, dtype=np.float64), source_path=path)
