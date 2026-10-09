"""LED array geometry: raster grid, spiral update order, misaligned wavevectors.

LEDs are indexed in raster order, row by row from the top-left corner, which
is also the order of the frames in the stored measurements. Positions are in
millimetres; ``x`` grows to the right and ``y`` upward.

A pose is the vector ``(x_mm, y_mm, theta_rad, height_mm)``: the global shift
of the array, its in-plane rotation and its distance to the sample. It is the
``zz`` vector of the MATLAB code, kept in the same order and units so the
annealing steps mean the same thing.
"""

from __future__ import annotations

import numpy as np


def led_grid(array_size: int, pitch_mm: float) -> np.ndarray:
    """Nominal ``(N, 2)`` positions ``(x, y)`` of a square array in raster order."""
    half = (array_size - 1) / 2
    row, col = np.divmod(np.arange(array_size ** 2), array_size)
    x = (col - half) * pitch_mm
    y = (half - row) * pitch_mm
    return np.stack([x, y], axis=1)


def spiral_order(array_size: int) -> np.ndarray:
    """Raster indices of the LEDs in the order they are used, spiralling out from the centre.

    A port of ``gseq.m``: the walk alternates between the row and the column
    direction with legs of length 1, 1, 2, 2, 3, 3, ..., so the first ``n**2``
    entries are exactly the central ``n x n`` block for every odd ``n``.
    """
    if array_size < 1 or array_size % 2 == 0:
        raise ValueError(f"array_size must be a positive odd number, got {array_size}")
    centre = (array_size + 1) // 2  # 1-based, as gseq(9, 9, 17)
    rows = [centre]
    cols = [centre]
    dx, dy = 1, -1
    step_x, step_y = 1, -1
    direction = 1
    counter = 0
    for _ in range(1, array_size ** 2):
        counter += 1
        if direction == 1:
            rows.append(rows[-1] + dx)
            cols.append(cols[-1])
            if counter == abs(step_x):
                counter = 0
                direction = -direction
                dx = -dx
                step_x = -step_x
                step_x = step_x + 1 if step_x > 0 else step_x - 1
        else:
            rows.append(rows[-1])
            cols.append(cols[-1] + dy)
            if counter == abs(step_y):
                counter = 0
                direction = -direction
                dy = -dy
                step_y = -step_y
                step_y = step_y + 1 if step_y > 0 else step_y - 1

    return (np.asarray(rows) - 1) * array_size + (np.asarray(cols) - 1)


def illumination_wavevectors(grid_mm: np.ndarray, pose, sample_offset_mm, k0: float):
    """Transverse wavevectors ``(kx, ky)`` of every LED, a port of ``fpmangle.m``.

    The array is rotated by ``theta`` about its own centre, shifted by
    ``(x_mm, y_mm)`` and placed ``height_mm`` below the sample. Each LED is seen
    from the centre of the reconstructed patch, ``sample_offset_mm = (x, y)``.
    """
    shift_x, shift_y, theta, height = (float(value) for value in pose)
    offset_x, offset_y = (float(value) for value in sample_offset_mm)
    x = grid_mm[:, 0]
    y = grid_mm[:, 1]

    x_led = np.cos(theta) * x + np.sin(theta) * y + shift_x
    y_led = -np.sin(theta) * x + np.cos(theta) * y + shift_y

    kx = k0 * np.sin(np.arctan((offset_x - x_led) / height))
    ky = k0 * np.sin(np.arctan((offset_y - y_led) / height))
    return kx, ky
