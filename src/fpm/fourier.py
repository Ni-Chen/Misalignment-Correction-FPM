"""Centred Fourier transforms and rounding half away from zero.

These are the primitives shared by the forward model and the solvers. They
operate on the trailing two (image) dimensions.
"""

from __future__ import annotations

import numpy as np
import torch


def FT(x, dim=(-2, -1)):
    """Centred 2-D Fourier transform, ``fftshift(fft2(x))``."""
    fft_x = torch.fft.fft2(x, dim=dim)
    x_prime = torch.fft.fftshift(fft_x, dim=dim)

    return x_prime


def iFT(x, dim=(-2, -1)):
    """Inverse of :func:`FT`, ``ifft2(ifftshift(x))``."""
    shift_x = torch.fft.ifftshift(x, dim=dim)
    fft_x = torch.fft.ifft2(shift_x, dim=dim)

    return fft_x


def round_half_away(x):
    """Round half away from zero.

    Python's ``round`` and ``np.round`` round half to even, so two centres a
    whole pixel apart could get windows zero or two pixels apart. With an
    even object size every nominal LED centre falls on a half pixel, so the
    windows need a rule that treats all of them alike.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)
