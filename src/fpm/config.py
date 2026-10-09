from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


METHOD_FILES = {
    "pie": "config_pie.yaml",
    "sa": "config_sa.yaml",
    "pcfpm": "config_pcfpm.yaml",
    "mcfpm": "config_mcfpm.yaml",
}

# ``data/`` holds acquisitions only. Everything a solver produces is written under ``notebooks/output/`` so a dataset directory stays exactly what was measured.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = PROJECT_ROOT / "data"
OUTPUT_ROOT = PROJECT_ROOT / "notebooks" / "output"


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    value = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(value, dict):
        raise TypeError(f"Configuration must be a mapping: {path}")
    return value


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def iter_dataset_dirs(data_root: str | Path):
    """Yield every dataset directory under ``data_root``.

    A dataset is any directory holding a ``config.yaml``. Related acquisitions
    are grouped one level deep -- ``data/root/{center,edge_upper_left,...}`` -- so the
    search recurses. Directories whose name starts with ``_`` (the template,
    shared assets) are not datasets.
    """
    root = Path(data_root)
    for path in sorted(root.rglob("config.yaml")):
        directory = path.parent
        relative = directory.relative_to(root)
        if any(part.startswith("_") for part in relative.parts):
            continue
        yield directory


def dataset_name(dataset_dir: str | Path, data_root: str | Path | None = None) -> str:
    """Name a dataset by its path relative to ``data_root``: ``root/center``, not ``center``."""
    path = Path(dataset_dir).expanduser().resolve()
    root = Path(data_root).expanduser().resolve() if data_root is not None else DATA_ROOT
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.name


def output_dir(
    dataset_dir: str | Path,
    method: str,
    *,
    data_root: str | Path | None = None,
    output_root: str | Path | None = None,
) -> Path:
    """Where one dataset/method result belongs: ``<output_root>/<name>/<method>``.

    The directory is not created here; ``ReconstructionResult.save`` does that.
    """
    root = (
        Path(output_root).expanduser().resolve()
        if output_root is not None
        else OUTPUT_ROOT
    )
    return root / dataset_name(dataset_dir, data_root) / str(method).lower()


def load_config(dataset_dir: str | Path, method: str) -> dict[str, Any]:
    """Load the dataset config and merge the selected method config.

    ``config.yaml`` describes the stored data and the optical system.
    ``config_<method>.yaml`` holds solver parameters only and is optional: a
    missing file, or a missing key, takes the default listed in
    ``data/_template``, which reproduces the MATLAB code. A method file may
    also override ``model.sample_offset_*`` to reproduce the geometry a
    particular MATLAB script used.
    """
    dataset_path = Path(dataset_dir).expanduser().resolve()
    method_name = str(method).lower()
    if method_name not in METHOD_FILES:
        choices = ", ".join(sorted(METHOD_FILES))
        raise ValueError(f"Unknown reconstruction method {method!r}; choose {choices}.")

    common_path = dataset_path / "config.yaml"
    method_path = dataset_path / METHOD_FILES[method_name]
    common = _read_yaml(common_path)
    specific = _read_yaml(method_path) if method_path.is_file() else {}

    model_override = specific.pop("model", None)
    merged = _deep_merge(common, {"reconstruction": specific})
    if model_override is not None:
        if not isinstance(model_override, dict):
            raise TypeError(f"'model' must be a mapping: {method_path}")
        merged = _deep_merge(merged, {"model": model_override})
    merged["method"] = method_name
    merged["dataset_dir"] = dataset_path
    merged["config_paths"] = {
        "common": common_path,
        "method": method_path if method_path.is_file() else None,
    }
    validate_config(merged)
    return merged


def validate_config(config: dict[str, Any]) -> None:
    data = config.get("data")
    if not isinstance(data, dict):
        raise KeyError("config is missing the 'data' mapping")
    for key in ("file", "key", "layout"):
        if key not in data:
            raise KeyError(f"config.data is missing '{key}'")
    layout = str(data["layout"]).upper()
    if len(layout) != 3 or set(layout) != {"I", "H", "W"}:
        raise ValueError("data.layout must be a permutation of 'IHW', for example IHW or HWI")

    for name in ("optics", "led"):
        if not isinstance(config.get(name), dict):
            raise KeyError(f"config is missing the '{name}' mapping")
    optics = config["optics"]
    led = config["led"]
    for key in ("wavelength_m", "na", "magnification", "pixel_size_sensor_m"):
        if key not in optics:
            raise KeyError(f"config.optics is missing '{key}'")
        if float(optics[key]) <= 0:
            raise ValueError(f"config.optics.{key} must be positive")
    for key in ("array_size", "pitch_m", "height_m"):
        if key not in led:
            raise KeyError(f"config.led is missing '{key}'")
    if int(led["array_size"]) % 2 == 0:
        raise ValueError("led.array_size must be odd: the spiral starts at the centre LED")

    reconstruction = config.get("reconstruction") or {}
    if int(reconstruction.get("upscale", 4)) < 1:
        raise ValueError("reconstruction.upscale must be at least 1")
