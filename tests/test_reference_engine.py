"""The ptychographic engine against an independent reference implementation.

``reference_sweeps`` is a plain NumPy version of the same sweep, written
without ``src/fpm``: it places every window by its own arithmetic, in 1-based
pixel coordinates counted from the first pixel of the spectrum, and turns them
into slices only when cutting. It shares nothing with the package but the LED
wavevectors and the update order, which ``test_geometry.py`` checks on their own.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from fpm import illumination_wavevectors, spiral_order
from support import cpu_model, simulated_dataset


def _round_half_away(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def reference_sweeps(stack, kx, ky, *, frame_size, upscale, na, k0, pixel_size, sweeps, round_centres):
    """``sweeps`` sweeps over every LED of ``stack`` (rows, columns, LEDs).

    With ``round_centres`` each centre is rounded to a whole pixel before its
    window is placed, as ``pie`` does.
    """
    size_hr = frame_size * upscale
    dk = 2 * np.pi / (pixel_size / upscale * size_hr)
    k_max = np.pi / pixel_size
    axis = np.linspace(-k_max, k_max, frame_size)
    kx_grid, ky_grid = np.meshgrid(axis, axis)
    support = ((kx_grid ** 2 + ky_grid ** 2) < (na * k0) ** 2).astype(np.float64)
    defocus = np.ones((frame_size, frame_size))
    array_size = int(round(np.sqrt(stack.shape[2])))
    order = spiral_order(array_size) + 1                    # 1-based LED numbers
    centre_led = (array_size ** 2 + 1) // 2
    pad = (size_hr - frame_size) // 2
    spectrum = np.pad(np.fft.fftshift(np.fft.fft2(stack[:, :, centre_led - 1])), pad)
    pupil = 1
    scale = (frame_size / size_hr) ** 2
    alpha = 1
    beta = 1e0
    errors = []
    for _ in range(sweeps):
        error = 0.0
        for led in order:
            centre_x = (size_hr + 1) / 2 + kx[led - 1] / dk
            centre_y = (size_hr + 1) / 2 + ky[led - 1] / dk
            if round_centres:
                centre_x = _round_half_away(centre_x)
                centre_y = _round_half_away(centre_y)
            first_row = int(_round_half_away(centre_y - (frame_size - 1) / 2))
            last_row = int(_round_half_away(centre_y + (frame_size - 1) / 2))
            first_col = int(_round_half_away(centre_x - (frame_size - 1) / 2))
            last_col = int(_round_half_away(centre_x + (frame_size - 1) / 2))
            rows = slice(first_row - 1, last_row)
            cols = slice(first_col - 1, last_col)
            intensity = stack[:, :, led - 1]

            predicted = scale * spectrum[rows, cols] * support * pupil * defocus
            field = np.fft.ifft2(np.fft.ifftshift(predicted))
            error += np.sum((scale * np.abs(field) - intensity) ** 2)
            field = (1 / scale) * np.sqrt(intensity) * np.exp(1j * np.angle(field))
            corrected = np.fft.fftshift(np.fft.fft2(field))
            spectrum[rows, cols] = spectrum[rows, cols] + (
                1 / np.max(np.abs(pupil)) * np.abs(pupil) * np.conj(pupil) * (corrected - predicted)
                / (np.abs(pupil) ** 2 + alpha)
            )
            patch = spectrum[rows, cols]
            pupil = pupil + (
                1 / np.max(np.abs(patch)) * np.abs(patch) * np.conj(patch) * (corrected - predicted)
                / (np.abs(patch) ** 2 + beta) * support
            )
        errors.append(error)
    return np.fft.ifft2(np.fft.ifftshift(spectrum)), pupil, errors


@pytest.mark.parametrize("round_centres", [False, True])
def test_engine_matches_the_reference(tmp_path, round_centres):
    model = cpu_model(simulated_dataset(tmp_path), "pie")
    pose = np.array([0.7, -0.4, np.deg2rad(0.5), 114.0])
    kx, ky = illumination_wavevectors(model.led_grid, pose, model.sample_offset_mm, model.k0)
    stack = np.transpose(model.intensity.numpy(), (1, 2, 0))

    expected_object, expected_pupil, expected_errors = reference_sweeps(
        stack, kx, ky, frame_size=model.size_lr, upscale=model.upscale, na=model.na, k0=model.k0,
        pixel_size=model.pixel_size_object, sweeps=2, round_centres=round_centres,
    )

    centres = model.spectrum_centres(pose)
    if round_centres:
        centres = np.sign(centres) * np.floor(np.abs(centres) + 0.5)
    spectrum, pupil, _, errors = model.reconstruct(centres, 2, cost="amplitude_vs_intensity")
    object_field = torch.fft.ifft2(torch.fft.ifftshift(spectrum)).numpy()

    scale = np.abs(expected_object).max()
    assert np.max(np.abs(object_field - expected_object)) / scale < 1e-10
    assert np.max(np.abs(pupil.numpy() - expected_pupil)) < 1e-10
    assert errors == pytest.approx(expected_errors, rel=1e-10)


def test_rounded_and_continuous_centres_cut_different_windows(tmp_path):
    """``pie`` rounds each centre before cutting its window; the misalignment solvers do not."""
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    centres = model.spectrum_centres()
    rounded = np.sign(centres) * np.floor(np.abs(centres) + 0.5)
    windows = np.array([model.window(c) for c in centres])
    rounded_windows = np.array([model.window(c) for c in rounded])
    assert np.any(windows != rounded_windows)
    # The aligned centre LED lands on (n - 1) / 2, a half pixel, which rounding moves up.
    assert tuple(rounded_windows[model.centre_frame] - windows[model.centre_frame]) == (1, 1)
