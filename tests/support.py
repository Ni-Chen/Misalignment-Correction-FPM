"""Shared by the test modules: small synthetic instruments and acquisitions.

The defaults are the real instrument with smaller frames (64 px instead of
128), which keeps a full 17 x 17 sweep under a tenth of a second while every
LED window still fits in the object spectrum.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from fpm import FPM
from fpm.config import METHOD_FILES, load_config
from fpm.simulation import load_object, simulate

ROOT = Path(__file__).resolve().parents[1]

OPTICS = {"wavelength_m": 6.297e-7, "na": 0.1, "magnification": 4.0, "pixel_size_sensor_m": 6.5e-6}
LED = {"array_size": 17, "pitch_m": 4.0e-3, "height_m": 113.5e-3}
DATA = {"file": "intensity.npy", "key": "intensity", "layout": "IHW"}

# A misaligned array: shifted by (1.5, -1.0) mm. With 64 px frames one spectrum
# pixel is about 0.69 mm of shift near the axis, so this is (2.2, -1.5) px.
TRUE_POSE = np.array([1.5, -1.0, 0.0, 113.5])


def write_yaml(path: Path, value: dict) -> None:
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def write_dataset(directory: Path, intensity, *, led: dict | None = None,
                  model: dict | None = None, **method_files: dict) -> Path:
    """Create a dataset directory holding ``intensity`` (frames, rows, columns) and return it.

    ``led`` extends the default array, ``model`` is written as given, and each
    method passed as a keyword (``mcfpm={...}`` and so on) gets its
    ``config_<method>.yaml``.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    np.save(directory / "intensity.npy", np.asarray(intensity, dtype=np.float64))
    common = {"data": dict(DATA), "optics": dict(OPTICS), "led": {**LED, **(led or {})}}
    if model is not None:
        common["model"] = model
    write_yaml(directory / "config.yaml", common)
    for method, settings in method_files.items():
        write_yaml(directory / METHOD_FILES[method], settings)
    return directory


def simulated_dataset(directory: Path, *, pose=TRUE_POSE, size_hr: int = 256, array_size: int = 17,
                      **method_files: dict) -> Path:
    """The test object seen through a (by default misaligned) array, written as a dataset."""
    placeholder = write_dataset(directory, np.zeros((array_size ** 2, 1, 1)), led={"array_size": array_size})
    config = load_config(placeholder, "pie")
    data = simulate(config, load_object(size_hr), pose, peak_intensity=3e4)
    return write_dataset(directory, data.values, led={"array_size": array_size}, **method_files)


def cpu_model(dataset_dir: Path, method: str) -> FPM:
    return FPM.from_directory(dataset_dir, method, device="cpu")
