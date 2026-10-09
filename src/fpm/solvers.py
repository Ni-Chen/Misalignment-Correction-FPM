"""Misalignment-correcting solvers built on the engine of :class:`~fpm.model.FPM`.

``pcfpm`` corrects every LED by simulated annealing and then pulls them back
onto a rigid array model; ``mcfpm``, the method of the paper, anneals the
array pose itself with the whole reconstruction as the cost. Both are
functions of a configured :class:`FPM` and are exposed as ``model.pcfpm()``
and ``model.mcfpm()``.

Default settings are the values in the MATLAB code, not those in the paper
where the two differ.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.optimize import least_squares

from .model import FPM, ReconstructionResult, _pose_for_humans
from .optimizers import simulated_annealing

# pie_solver_coarse.m divides its error by this before handing it to the
# annealer. The acceptance of an uphill step depends on cost / temperature,
# so the constant is part of the search, not cosmetic.
COARSE_COST_SCALE = 1e4

PCFPM_BOUNDS = {"x_mm": [-4.0, 4.0], "y_mm": [-4.0, 4.0], "rotation_deg": [-2.0, 2.0], "height_mm": [0.0, 6.5]}
MCFPM_BOUNDS = {"x_mm": [-4.0, 4.0], "y_mm": [-4.0, 4.0], "rotation_deg": [0.0, 0.0], "height_mm": [0.0, 0.0]}


def fit_led_pose(model: FPM, centres: np.ndarray, pose0, lower, upper, frames) -> np.ndarray:
    """Least-squares array pose whose spectrum centres best match ``centres`` on ``frames``.

    The regression step of pcFPM (``fpm_fit.m`` under ``fmincon``). A
    dimension whose bounds coincide is held at that value.
    """
    pose0 = np.clip(np.asarray(pose0, dtype=np.float64), lower, upper)
    free = lower < upper
    frames = np.asarray(frames)

    def residuals(free_values):
        pose = pose0.copy()
        pose[free] = free_values
        return (model.spectrum_centres(pose)[frames] - centres[frames]).ravel()

    if not free.any():
        return pose0
    fit = least_squares(residuals, pose0[free], bounds=(lower[free], upper[free]), x_scale="jac")
    pose = pose0.copy()
    pose[free] = fit.x
    return pose


def pcfpm(model: FPM) -> ReconstructionResult:
    """Position-correcting FPM, the comparison method of ``pie_solver_led.m``.

    Every outer iteration restarts the object and pupil and runs
    ``sweeps_per_iteration`` sweeps; the last of them anneals each LED's
    centre within the search radius before using it. After every sweep the
    array pose is refitted to the current centres, and the centres are
    replaced by those of the fitted array. While ``outer < fit_warmup_iterations``
    only the central ``fit_array_size x fit_array_size`` LEDs enter the fit.

    The search radius is ``search_radius_px / search_radius_divisor`` for the
    first ``divided_search_iterations`` outer iterations and
    ``search_radius_px`` afterwards, exactly as the MATLAB code resets and
    divides ``dkshift``.

    ``result.loss`` holds the intensity misfit of every sweep and
    ``result.led_pose_history`` the fitted pose after every outer iteration.
    """
    recon = model.recon
    iterations = int(recon.get("iterations", 9))
    sweeps = int(recon.get("sweeps_per_iteration", 11))
    radius = float(recon.get("search_radius_px", 12.0))
    divisor = float(recon.get("search_radius_divisor", 1.3))
    divided_iterations = int(recon.get("divided_search_iterations", 8))
    fit_size = int(recon.get("fit_array_size", 5))
    warmup = int(recon.get("fit_warmup_iterations", 5))
    lower, upper = model.pose_bounds(recon.get("pose_bounds", PCFPM_BOUNDS), model.nominal_pose)
    rng = np.random.default_rng(int(recon["seed"]))

    pose = model.nominal_pose
    centres = model.spectrum_centres(pose)
    errors: list[float] = []
    pose_history: list[list[float]] = []
    for outer in range(iterations):
        spectrum = model.initial_spectrum()
        pupil = model.initial_pupil()
        search_radius = radius / divisor if outer < divided_iterations else radius
        fit_frames = model.scan_order[: fit_size ** 2] if outer < warmup else model.scan_order
        for inner in range(sweeps):
            searching = inner == sweeps - 1
            pupil, error = model.sweep(spectrum, pupil, centres, model.scan_order,
                                       search_radius=search_radius if searching else None, rng=rng)
            errors.append(error)
            pose = fit_led_pose(model, centres, pose, lower, upper, fit_frames)
            centres = model.spectrum_centres(pose)
        pose_history.append(_pose_for_humans(pose).tolist())
        if int(recon["print_every"]):
            print(f"pcfpm {outer + 1}/{iterations}: pose {np.round(pose_history[-1], 4)}", flush=True)

    return model._result("pcfpm", spectrum, pupil, centres, errors, pose=pose,
                         led_pose_history=pose_history)


def mcfpm(model: FPM) -> ReconstructionResult:
    """Misalignment-correcting FPM (``main.m`` with ``pie_solver_coarse.m`` and ``pie_solver_fine.m``).

    *Coarse.* The array pose ``(x, y, theta, height)`` is annealed within
    ``pose_bounds``; each candidate is scored by a reconstruction of
    ``coarse_iterations`` sweeps from scratch over the central
    ``coarse_array_size x coarse_array_size`` LEDs, whose last-sweep error,
    divided by :data:`COARSE_COST_SCALE`, is the cost. ``coarse_cost``
    chooses that error: ``legacy`` compares the modelled amplitude with the
    measured intensity, as the MATLAB code does, and ``intensity`` compares
    intensities, as Eq. 15 of the paper states.

    *Fine.* A reconstruction of ``fine_iterations`` sweeps over all LEDs with
    the annealed pose. From sweep ``fine_plain_sweeps`` on, each LED's centre
    is additionally annealed within ``fine_search_radius_px`` -- the paper's
    mcFPM III. With the MATLAB setting of two sweeps that never happens.

    The defaults are the MATLAB ones: all 17 x 17 LEDs in the cost (the
    paper's mcFPM II; ``coarse_array_size: 5`` is mcFPM I), 10 coarse sweeps
    (the paper says 5), rotation and height held at their nominal values.

    ``result.anneal_history`` holds the best coarse cost after every
    annealing iteration and ``result.loss`` the intensity misfit of every fine sweep.
    """
    recon = model.recon
    coarse_iterations = int(recon.get("coarse_iterations", 10))
    coarse_size = int(recon.get("coarse_array_size", model.array_size))
    coarse_cost = str(recon.get("coarse_cost", "legacy"))
    fine_iterations = int(recon.get("fine_iterations", 2))
    fine_plain_sweeps = int(recon.get("fine_plain_sweeps", 9))
    fine_radius = float(recon.get("fine_search_radius_px", 3.0))
    print_every = int(recon["print_every"])
    lower, upper = model.pose_bounds(recon.get("pose_bounds", MCFPM_BOUNDS), model.nominal_pose)
    coarse_frames = model.scan_order[: coarse_size ** 2]

    def cost(pose) -> float:
        _, _, _, errors = model.reconstruct(model.spectrum_centres(pose), coarse_iterations,
                                            order=coarse_frames, cost=coarse_cost)
        return errors[-1] / COARSE_COST_SCALE

    started = time.perf_counter()

    def report(iteration, best_pose, best_cost):
        if print_every and iteration % print_every == 0:
            print(f"mcfpm anneal {iteration}: cost {best_cost:.6g}, "
                  f"pose {np.round(_pose_for_humans(best_pose), 4)} "
                  f"({time.perf_counter() - started:.0f} s)", flush=True)

    anneal = simulated_annealing(
        cost,
        model.nominal_pose,
        lower,
        upper,
        max_iterations=int(recon.get("anneal_max_iterations", 100)),
        function_tolerance=float(recon.get("anneal_function_tolerance", 1e-3)),
        initial_temperature=float(recon.get("anneal_initial_temperature", 100.0)),
        seed=int(recon["seed"]),
        callback=report,
    )

    rng = np.random.default_rng(int(recon["seed"]) + 1)
    spectrum, pupil, centres, errors = model.reconstruct(
        model.spectrum_centres(anneal.x), fine_iterations,
        plain_sweeps=fine_plain_sweeps, search_radius=fine_radius, rng=rng,
    )
    result = model._result("mcfpm", spectrum, pupil, centres, errors, pose=anneal.x,
                           anneal_history=list(anneal.history))
    result.metadata["anneal_iterations"] = anneal.iterations
    result.metadata["anneal_function_evaluations"] = anneal.function_evaluations
    result.metadata["seconds"] = time.perf_counter() - started
    return result
