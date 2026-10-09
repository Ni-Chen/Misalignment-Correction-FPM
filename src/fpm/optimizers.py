"""Simulated annealing in a box.

Every misalignment search -- the global LED pose of mcFPM and the per-LED
spectrum positions of SA and pcFPM -- runs on this annealer:

* a step of length ``T`` in a uniformly random direction;
* a step leaving the box is projected onto it and then replaced by a random
  convex combination of that projection and the current point, so a
  dimension with equal bounds stays fixed;
* a worse point is accepted with probability ``1 / (1 + exp(delta / T))``;
* ``T = T0 * 0.95**k`` with ``T0 = 100``;
* stopping after ``max_iterations``, ``3000 * n`` function evaluations, or
  once the best value has changed by less than ``function_tolerance`` per
  iteration on average over the last ``500 * n`` iterations.
"""

from __future__ import annotations

from typing import Callable, NamedTuple

import numpy as np
from scipy.special import expit


class AnnealingResult(NamedTuple):
    x: np.ndarray               # best point found
    fun: float                  # its cost
    iterations: int
    function_evaluations: int
    history: list[float]        # best cost before the first and after every iteration


def simulated_annealing(
    cost: Callable[[np.ndarray], float],
    x0,
    lower,
    upper,
    *,
    max_iterations: int | None = None,
    max_function_evaluations: int | None = None,
    function_tolerance: float = 1e-6,
    stall_iterations: int | None = None,
    initial_temperature: float = 100.0,
    seed: int | None = None,
    callback: Callable[[int, np.ndarray, float], None] | None = None,
) -> AnnealingResult:
    """Minimise ``cost`` over the box ``[lower, upper]`` starting from ``x0``.

    ``callback(iteration, best_x, best_cost)`` is called after every iteration.
    """
    x = np.array(x0, dtype=np.float64).ravel()
    lower = np.broadcast_to(np.asarray(lower, dtype=np.float64), x.shape)
    upper = np.broadcast_to(np.asarray(upper, dtype=np.float64), x.shape)
    if np.any(lower > upper):
        raise ValueError("every lower bound must not exceed its upper bound")
    if np.any((x < lower) | (x > upper)):
        raise ValueError("x0 lies outside the bounds")

    n = x.size
    max_iterations = np.inf if max_iterations is None else max_iterations
    max_function_evaluations = 3000 * n if max_function_evaluations is None else max_function_evaluations
    stall_iterations = 500 * n if stall_iterations is None else stall_iterations
    rng = np.random.default_rng(seed)

    current_cost = float(cost(x))
    evaluations = 1
    best_x, best_cost = x.copy(), current_cost
    history = [best_cost]
    temperature = initial_temperature
    iteration = 0
    tiny = np.finfo(np.float64).tiny

    while iteration < max_iterations and evaluations < max_function_evaluations:
        if iteration >= stall_iterations:
            average_change = abs(history[-1 - stall_iterations] - history[-1]) / stall_iterations
            if average_change < function_tolerance:
                break

        direction = rng.standard_normal(n)
        direction /= np.linalg.norm(direction)
        candidate = x + temperature * direction
        if np.any((candidate < lower) | (candidate > upper)):
            weight = rng.random()
            candidate = weight * np.clip(candidate, lower, upper) + (1 - weight) * x

        candidate_cost = float(cost(candidate))
        evaluations += 1
        delta = candidate_cost - current_cost
        if delta < 0 or expit(-delta / max(temperature, tiny)) > rng.random():
            x, current_cost = candidate, candidate_cost
        if candidate_cost < best_cost:
            best_x, best_cost = candidate.copy(), candidate_cost

        iteration += 1
        temperature = initial_temperature * 0.95 ** iteration
        history.append(best_cost)
        if callback is not None:
            callback(iteration, best_x, best_cost)

    return AnnealingResult(best_x, best_cost, iteration, evaluations, history)
