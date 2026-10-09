"""The ptychographic engine against a line-by-line transcription of the MATLAB loop.

No MATLAB output could be produced for this port, so the reference here is
``pie_solver_fine.m`` (and the rounding variant of ``pie_solver.m``) written
out in NumPy with MATLAB's 1-based indices, column-major frame stack and
expression order. It shares nothing with ``src/fpm`` except the LED
wavevectors, which ``test_geometry.py`` checks against ``fpmangle.m`` on its
own.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from fpm import illumination_wavevectors, spiral_order
from support import cpu_model, simulated_dataset


def _mround(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def matlab_sweeps(imSeqLowRes, kx, ky, *, m1, mag_image, NA, k0, spsize, loop, round_centres):
    """pie_solver_fine.m with loop < 10 (round_centres=False) or pie_solver.m (True)."""
    n1 = m1
    m = m1 * mag_image
    n = m
    psize = spsize / mag_image
    dkx = 2 * np.pi / (psize * n)
    dky = 2 * np.pi / (psize * m)
    kmax = np.pi / spsize
    axis = np.linspace(-kmax, kmax, n1)
    kxm, kym = np.meshgrid(axis, axis)
    CTF = ((kxm ** 2 + kym ** 2) < (NA * k0) ** 2).astype(np.float64)
    H = np.ones((m1, n1))
    arraysize = int(round(np.sqrt(imSeqLowRes.shape[2])))
    seq = spiral_order(arraysize) + 1                       # gseq(9, 9, arraysize)
    centre = (arraysize ** 2 + 1) // 2                      # 145
    pad = (m - m1) // 2
    objectRecoverFT = np.pad(np.fft.fftshift(np.fft.fft2(imSeqLowRes[:, :, centre - 1])), pad)
    pupil = 1
    errorpie = []
    alpha = 1
    beta = 1e0
    for tt in range(1, loop + 1):
        error = 0.0
        for i3 in range(1, arraysize ** 2 + 1):
            i2 = seq[i3 - 1]
            kxc = (n + 1) / 2 + kx[i2 - 1] / dkx
            kyc = (n + 1) / 2 + ky[i2 - 1] / dky
            if round_centres:
                kxc = _mround(kxc)
                kyc = _mround(kyc)
            kyl = int(_mround(kyc - (m1 - 1) / 2))
            kyh = int(_mround(kyc + (m1 - 1) / 2))
            kxl = int(_mround(kxc - (n1 - 1) / 2))
            kxh = int(_mround(kxc + (n1 - 1) / 2))
            rows = slice(kyl - 1, kyh)
            cols = slice(kxl - 1, kxh)
            lowResFT_1 = (m1 / m) ** 2 * objectRecoverFT[rows, cols] * CTF * pupil * H
            im_lowRes = np.fft.ifft2(np.fft.ifftshift(lowResFT_1))
            error += np.sum(((m1 / m) ** 2 * np.abs(im_lowRes) - imSeqLowRes[:, :, i2 - 1]) ** 2)
            im_lowRes = (m / m1) ** 2 * np.sqrt(imSeqLowRes[:, :, i2 - 1]) * np.exp(1j * np.angle(im_lowRes))
            lowResFT_2 = np.fft.fftshift(np.fft.fft2(im_lowRes))
            objectRecoverFT[rows, cols] = objectRecoverFT[rows, cols] + (
                1 / np.max(np.abs(pupil)) * np.abs(pupil) * np.conj(pupil) * (lowResFT_2 - lowResFT_1)
                / (np.abs(pupil) ** 2 + alpha)
            )
            patch = objectRecoverFT[rows, cols]
            pupil = pupil + (
                1 / np.max(np.abs(patch)) * np.abs(patch) * np.conj(patch) * (lowResFT_2 - lowResFT_1)
                / (np.abs(patch) ** 2 + beta) * CTF
            )
        errorpie.append(error)
    return np.fft.ifft2(np.fft.ifftshift(objectRecoverFT)), pupil, errorpie


@pytest.mark.parametrize("round_centres", [False, True])
def test_engine_matches_the_matlab_loop(tmp_path, round_centres):
    model = cpu_model(simulated_dataset(tmp_path), "pie")
    pose = np.array([0.7, -0.4, np.deg2rad(0.5), 114.0])
    kx, ky = illumination_wavevectors(model.led_grid, pose, model.sample_offset_mm, model.k0)
    stack = np.transpose(model.intensity.numpy(), (1, 2, 0))          # HWI, as in the .mat files

    expected_object, expected_pupil, expected_errors = matlab_sweeps(
        stack, kx, ky, m1=model.size_lr, mag_image=model.upscale, NA=model.na, k0=model.k0,
        spsize=model.pixel_size_object, loop=2, round_centres=round_centres,
    )

    centres = model.spectrum_centres(pose)
    if round_centres:
        centres = np.sign(centres) * np.floor(np.abs(centres) + 0.5)
    spectrum, pupil, _, errors = model.reconstruct(centres, 2, cost="legacy")
    object_field = torch.fft.ifft2(torch.fft.ifftshift(spectrum)).numpy()

    scale = np.abs(expected_object).max()
    assert np.max(np.abs(object_field - expected_object)) / scale < 1e-10
    assert np.max(np.abs(pupil.numpy() - expected_pupil)) < 1e-10
    assert errors == pytest.approx(expected_errors, rel=1e-10)


def test_rounded_and_continuous_centres_cut_different_windows(tmp_path):
    """pie_solver.m rounds each centre before cutting its window; the misalignment solvers do not."""
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    centres = model.spectrum_centres()
    rounded = np.sign(centres) * np.floor(np.abs(centres) + 0.5)
    windows = np.array([model.window(c) for c in centres])
    rounded_windows = np.array([model.window(c) for c in rounded])
    assert np.any(windows != rounded_windows)
    # The aligned centre LED lands on (n - 1) / 2, a half pixel, which rounding moves up.
    assert tuple(rounded_windows[model.centre_frame] - windows[model.centre_frame]) == (1, 1)
