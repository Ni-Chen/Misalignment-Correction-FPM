"""The simulated annealer and the pcFPM array fit."""

from __future__ import annotations

import numpy as np
import pytest

from fpm import fit_led_pose, simulated_annealing
from support import cpu_model, simulated_dataset


def bowl(x):
    return float(np.sum((np.asarray(x) - np.array([1.2, -0.7])) ** 2))


def test_annealing_finds_the_minimum_of_a_bowl():
    result = simulated_annealing(bowl, [0.0, 0.0], [-4.0, -4.0], [4.0, 4.0],
                                 max_iterations=400, seed=0)
    assert result.x == pytest.approx([1.2, -0.7], abs=0.05)
    assert result.iterations == 400
    assert result.function_evaluations == 401
    assert len(result.history) == 401
    assert np.all(np.diff(result.history) <= 0)


def test_annealing_is_reproducible_and_stays_in_bounds():
    visited = []

    def cost(x):
        visited.append(np.array(x))
        return bowl(x)

    first = simulated_annealing(cost, [0.0, 0.0], [-1.0, -1.0], [1.0, 1.0], max_iterations=50, seed=3)
    second = simulated_annealing(bowl, [0.0, 0.0], [-1.0, -1.0], [1.0, 1.0], max_iterations=50, seed=3)
    assert np.array_equal(first.x, second.x)
    visited = np.array(visited)
    assert visited.min() >= -1.0 and visited.max() <= 1.0
    assert first.fun < bowl([0.0, 0.0])


def test_equal_bounds_hold_a_parameter_fixed():
    visited = []

    def cost(x):
        visited.append(np.array(x))
        return bowl(x[:2])

    simulated_annealing(cost, [0.0, 0.0, 113.5], [-4.0, -4.0, 113.5], [4.0, 4.0, 113.5],
                        max_iterations=30, seed=0)
    assert np.all(np.array(visited)[:, 2] == 113.5)


def test_stall_criterion_stops_a_flat_search():
    result = simulated_annealing(lambda x: 1.0, [0.0], [-1.0], [1.0],
                                 function_tolerance=1e-3, stall_iterations=20, seed=0)
    assert result.iterations == 20


def test_array_fit_recovers_the_pose_of_exact_centres(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pcfpm")
    truth = np.array([0.8, -1.3, np.deg2rad(1.2), 116.0])
    lower = model.nominal_pose + np.array([-4, -4, np.deg2rad(-2), 0])
    upper = model.nominal_pose + np.array([4, 4, np.deg2rad(2), 6.5])
    pose = fit_led_pose(model, model.spectrum_centres(truth), model.nominal_pose, lower, upper,
                        model.scan_order)
    assert pose == pytest.approx(truth, abs=1e-6)
