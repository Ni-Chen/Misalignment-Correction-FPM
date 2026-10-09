"""Simulated acquisitions with a known object and a known LED misalignment.

The frames are made with the reconstruction's own forward model (see the
module docstring of :mod:`fpm.model`), so a solver that finds the right pose
can fit them exactly. This is the reference the solvers are tested against.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch
from scipy import ndimage

from .config import DATA_ROOT
from .data import IntensityData
from .fourier import FT, iFT
from .model import FPM

ASSETS = DATA_ROOT / "_assets"


def _read_gray(path: Path) -> np.ndarray:
    import tifffile

    image = np.asarray(tifffile.imread(path), dtype=np.float64)
    if image.ndim == 3:
        # ITU-R BT.601 luma weights
        image = image[..., :3] @ np.array([0.2989, 0.5870, 0.1140])
    return image


def _resize(image: np.ndarray, size: int) -> np.ndarray:
    return ndimage.zoom(image, (size / image.shape[0], size / image.shape[1]), order=1)


def load_object(size: int, *, amplitude_path: str | Path = ASSETS / "4.2.03.tiff",
                phase_path: str | Path = ASSETS / "5.2.09.tiff",
                phase_range: float = 0.5 * np.pi) -> np.ndarray:
    """A complex test object built from two standard images.

    The amplitude is the square root of the grey-scale mandrill, the phase
    the aerial image stretched to ``[0, phase_range]``, both resized to
    ``size x size``.
    """
    amplitude = np.sqrt(_resize(_read_gray(Path(amplitude_path)), size).clip(min=0))
    phase = _resize(_read_gray(Path(phase_path)), size)
    phase = (phase - phase.min()) / (phase.max() - phase.min()) * phase_range
    return amplitude * np.exp(1j * phase)


def simulate(config: dict[str, Any], object_field: np.ndarray, pose=None, *,
             peak_intensity: float | None = None) -> IntensityData:
    """Frames the instrument of ``config`` records of ``object_field`` with the array at ``pose``.

    ``config`` is a merged configuration as :func:`fpm.config.load_config`
    returns; its ``reconstruction.upscale`` sets the frame size as
    ``object size / upscale``. ``pose`` is ``(x_mm, y_mm, theta_rad, height_mm)``
    and defaults to the nominal array. ``peak_intensity`` rescales the stack
    so that its brightest pixel has that value, which matters for the
    ``amplitude_vs_intensity`` mcFPM cost: it mixes amplitude and intensity and therefore
    depends on the absolute level of the data.
    """
    object_field = np.asarray(object_field)
    size_hr = object_field.shape[-1]
    upscale = int((config.get("reconstruction") or {}).get("upscale", FPM.DEFAULTS["upscale"]))
    if size_hr % upscale:
        raise ValueError(f"object size {size_hr} is not a multiple of upscale {upscale}")
    size_lr = size_hr // upscale
    n_frames = int(config["led"]["array_size"]) ** 2
    placeholder = IntensityData(np.zeros((n_frames, size_lr, size_lr)), None)
    model = FPM(config, placeholder, device="cpu")

    spectrum = FT(torch.as_tensor(object_field, dtype=model.complex_dtype))
    centres = model.spectrum_centres(pose)
    frames = np.empty((n_frames, size_lr, size_lr))
    for frame, centre in enumerate(centres):
        y0, x0 = model.window(centre)
        patch = spectrum[y0:y0 + size_lr, x0:x0 + size_lr]
        field = iFT(model._scale * patch * model.transfer_function * model.defocus_kernel)
        frames[frame] = (model._scale * field.abs()).square().numpy()

    if peak_intensity is not None:
        frames *= peak_intensity / frames.max()
    return IntensityData(values=frames, source_path=None)
