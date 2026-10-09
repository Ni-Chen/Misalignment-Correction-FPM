"""LED geometry, update order and rounding: the ports of gseq.m, fpmangle.m and MATLAB's round."""

from __future__ import annotations

import numpy as np
import pytest

from fpm import illumination_wavevectors, led_grid, matlab_round, spiral_order
from support import cpu_model, simulated_dataset


def test_matlab_round_rounds_half_away_from_zero():
    values = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 192.5, 1.4999])
    assert matlab_round(values).tolist() == [-3, -2, -1, 1, 2, 3, 193, 1]


@pytest.mark.parametrize("array_size", [1, 3, 5, 17])
def test_spiral_visits_every_led_once_with_unit_steps(array_size):
    order = spiral_order(array_size)
    assert sorted(order.tolist()) == list(range(array_size ** 2))
    rows, cols = np.divmod(order, array_size)
    assert np.all(np.abs(np.diff(rows)) + np.abs(np.diff(cols)) == 1)


def test_spiral_matches_gseq_and_covers_central_blocks():
    order = spiral_order(17)
    # gseq(9, 9, 17) in MATLAB, minus one
    assert order[:9].tolist() == [144, 161, 160, 143, 126, 127, 128, 145, 162]
    for size in (1, 3, 5, 7):
        rows, cols = np.divmod(order[: size ** 2], 17)
        assert rows.min() == cols.min() == 8 - size // 2
        assert rows.max() == cols.max() == 8 + size // 2


def test_spiral_rejects_even_arrays():
    with pytest.raises(ValueError):
        spiral_order(4)


def test_led_grid_is_raster_from_top_left():
    grid = led_grid(3, 4.0)
    assert grid[0].tolist() == [-4.0, 4.0]
    assert grid[2].tolist() == [4.0, 4.0]
    assert grid[4].tolist() == [0.0, 0.0]
    assert grid[8].tolist() == [4.0, -4.0]


def test_wavevectors_follow_fpmangle():
    grid = led_grid(17, 4.0)
    k0 = 2 * np.pi / 6.297e-7
    pose = (0.3, -0.2, np.deg2rad(1.0), 113.5)
    offset = (0.5, -0.7)
    kx, ky = illumination_wavevectors(grid, pose, offset, k0)
    # One LED by hand: raster index 8 * 17 + 12 is (x, y) = (16, 0) mm.
    x, y = 16.0, 0.0
    theta = np.deg2rad(1.0)
    x_led = np.cos(theta) * x + np.sin(theta) * y + 0.3
    y_led = -np.sin(theta) * x + np.cos(theta) * y - 0.2
    index = 8 * 17 + 12
    assert kx[index] == pytest.approx(k0 * np.sin(np.arctan((0.5 - x_led) / 113.5)), rel=1e-14)
    assert ky[index] == pytest.approx(k0 * np.sin(np.arctan((-0.7 - y_led) / 113.5)), rel=1e-14)


def test_centre_led_of_an_aligned_array_sits_between_the_middle_pixels(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    centres = model.spectrum_centres()
    origin = (model.size_hr - 1) / 2
    assert centres[model.centre_frame].tolist() == [origin, origin]
    # (n + 1) / 2 - (m1 - 1) / 2 in MATLAB's 1-based indices, i.e. the block whose
    # middle pixel is the zero frequency of fftshift.
    assert model.window(centres[model.centre_frame]) == ((model.size_hr - model.size_lr) // 2,) * 2
    # Mirror-symmetric array, mirror-symmetric centres.
    assert np.allclose(centres[::-1] - origin, -(centres - origin))


def test_window_refuses_to_leave_the_spectrum(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    with pytest.raises(ValueError, match="outside"):
        model.window((model.size_hr, model.size_hr / 2))
