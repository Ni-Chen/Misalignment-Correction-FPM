"""Every solver on a simulated acquisition whose misalignment is known."""

from __future__ import annotations

import numpy as np
import pytest

from fpm.simulation import load_object
from support import TRUE_POSE, cpu_model, simulated_dataset

# Fast settings: mcFPM I (central 5 x 5 LEDs) with the paper's J = 5 sweeps.
MCFPM_FAST = {"coarse_array_size": 5, "coarse_iterations": 5, "print_every": 0}
SEARCH_FAST = {"search_stall_iterations": 15, "search_max_function_evaluations": 40, "print_every": 0}


def _correlation(estimate: np.ndarray, truth: np.ndarray) -> float:
    estimate = estimate - estimate.mean()
    truth = truth - truth.mean()
    return float(np.sum(estimate * truth) / np.sqrt(np.sum(estimate ** 2) * np.sum(truth ** 2)))


def test_pie_reconstructs_an_aligned_acquisition(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    result = model.run()
    assert result.loss[-1] < result.loss[0] / 10
    # The truth has detail beyond the synthetic NA (0.37 of the 0.78 the object grid holds).
    truth = np.abs(load_object(model.size_hr))
    assert _correlation(np.abs(result.object.numpy()), truth) > 0.7
    assert result.led_pose.tolist() == [0.0, 0.0, 0.0, 113.5]


def test_rounding_the_centres_costs_pie_accuracy(tmp_path):
    """pie_solver.m rounds the centres; on an aligned acquisition that alone raises the misfit.

    The nominal centre of an even-sized spectrum is a half pixel, and rounding
    it moves the window one pixel off the zero frequency -- a half-pixel bias
    the conventional-FPM baseline of the MATLAB comparisons carries. ``pie``
    keeps it to reproduce that baseline.
    """
    model = cpu_model(simulated_dataset(tmp_path, pose=None), "pie")
    rounded = model.run().loss[-1]
    continuous = model.reconstruct(model.spectrum_centres(), 9)[3][-1]
    assert rounded > 1.3 * continuous


def test_misalignment_degrades_pie(tmp_path):
    aligned = cpu_model(simulated_dataset(tmp_path / "aligned", pose=None), "pie").run()
    shifted = cpu_model(simulated_dataset(tmp_path / "shifted"), "pie").run()
    assert shifted.loss[-1] > 3 * aligned.loss[-1]


@pytest.mark.parametrize("cost", ["legacy", "intensity"])
def test_mcfpm_recovers_the_array_shift(tmp_path, cost):
    """The annealed pose cuts every LED's window within a pixel of where the truth cuts it.

    Windows sit on whole pixels, so the cost is flat across the poses that cut
    the same windows and the pose itself is only determined to about a
    pixel's worth of shift (0.69 mm here).
    """
    dataset = simulated_dataset(tmp_path, mcfpm={**MCFPM_FAST, "coarse_cost": cost})
    model = cpu_model(dataset, "mcfpm")
    result = model.run()
    pose = np.array([*result.led_pose[:2], 0.0, 113.5])
    found = np.array([model.window(c) for c in model.spectrum_centres(pose)])
    truth = np.array([model.window(c) for c in model.spectrum_centres(TRUE_POSE)])
    nominal = np.array([model.window(c) for c in model.spectrum_centres()])
    assert np.abs(found - truth).max() <= 1
    assert np.abs(nominal - truth).max() >= 2
    assert result.led_pose[2:].tolist() == [0.0, 113.5]
    assert len(result.anneal_history) == 101
    assert result.anneal_history[-1] < result.anneal_history[0]


def test_intensity_cost_separates_poses_more_sharply_than_legacy(tmp_path):
    """The MATLAB cost compares an amplitude with an intensity; Eq. 15 compares intensities."""
    model = cpu_model(simulated_dataset(tmp_path), "mcfpm")
    order = model.scan_order[:25]

    def cost(pose, kind):
        return model.reconstruct(model.spectrum_centres(pose), 5, order=order, cost=kind)[3][-1]

    for kind, minimum_ratio in (("legacy", 1.0002), ("intensity", 2.0)):
        assert cost(model.nominal_pose, kind) / cost(TRUE_POSE, kind) > minimum_ratio


def test_mcfpm_iii_searches_each_led_after_the_plain_sweeps(tmp_path):
    dataset = simulated_dataset(tmp_path, mcfpm={
        **MCFPM_FAST, **SEARCH_FAST, "anneal_max_iterations": 5,
        "fine_iterations": 2, "fine_plain_sweeps": 1,
    })
    result = cpu_model(dataset, "mcfpm").run()
    model = cpu_model(dataset, "mcfpm")
    annealed = model.spectrum_centres(np.array([*result.led_pose[:2], 0.0, 113.5]))
    moved = np.linalg.norm(result.spectrum_centres - annealed, axis=1)
    assert moved.max() > 0
    assert moved.max() <= np.sqrt(2) * 3.0 + 1e-9


def test_sa_moves_leds_independently_within_the_radius(tmp_path):
    dataset = simulated_dataset(tmp_path, sa={**SEARCH_FAST, "iterations": 2, "plain_sweeps": 1})
    result = cpu_model(dataset, "sa").run()
    start = np.sign(result.initial_spectrum_centres) * np.floor(np.abs(result.initial_spectrum_centres) + 0.5)
    moved = np.abs(result.spectrum_centres - start)
    assert moved.max() > 0
    assert moved.max() <= 10.0
    assert result.led_pose is None


def test_pcfpm_keeps_the_leds_on_a_rigid_array(tmp_path):
    dataset = simulated_dataset(tmp_path, pcfpm={
        **SEARCH_FAST, "iterations": 2, "sweeps_per_iteration": 2,
        "divided_search_iterations": 1, "fit_warmup_iterations": 1,
    })
    model = cpu_model(dataset, "pcfpm")
    result = model.run()
    pose = np.array(result.led_pose, dtype=np.float64)
    pose[2] = np.deg2rad(pose[2])
    assert result.spectrum_centres == pytest.approx(model.spectrum_centres(pose))
    assert len(result.led_pose_history) == 2
    assert len(result.loss) == 4
    lower = np.array([-4, -4, -2, 113.5])
    upper = np.array([4, 4, 2, 120.0])
    assert np.all(result.led_pose >= lower - 1e-9) and np.all(result.led_pose <= upper + 1e-9)
