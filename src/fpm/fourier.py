"""Centred Fourier transforms and MATLAB-compatible rounding.

These are the primitives shared by the forward model and the solvers. They
operate on the trailing two (image) dimensions.
"""

from __future__ import annotations

import numpy as np
import torch


def FT(x, dim=(-2, -1)):
    """Centred 2-D Fourier transform, MATLAB's ``fftshift(fft2(x))``."""
    fft_x = torch.fft.fft2(x, dim=dim)
    x_prime = torch.fft.fftshift(fft_x, dim=dim)

    return x_prime


def iFT(x, dim=(-2, -1)):
    """Inverse of :func:`FT`, MATLAB's ``ifft2(ifftshift(x))``."""
    shift_x = torch.fft.ifftshift(x, dim=dim)
    fft_x = torch.fft.ifft2(shift_x, dim=dim)

    return fft_x


def matlab_round(x):
    """Round half away from zero, as MATLAB's ``round`` does.

    Python's ``round`` and ``np.round`` round half to even, which moves a
    spectrum window by one pixel whenever a centre falls exactly on a half
    pixel -- and with an even object size the nominal centre of every LED
    window does.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)
