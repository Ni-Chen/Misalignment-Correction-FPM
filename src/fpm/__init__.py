"""Config-driven Fourier ptychography with LED misalignment correction.

The public surface is deliberately small: point :class:`FPM` at a dataset
directory, pick a method, call ``run()``.

    from fpm import FPM

    model = FPM.from_directory("data/usaf/edge_upper_right", method="mcfpm")
    result = model.run()
"""

from .config import DATA_ROOT, OUTPUT_ROOT, dataset_name, load_config, output_dir
from .data import IntensityData, load_intensity
from .fourier import FT, iFT, matlab_round
from .illumination import illumination_wavevectors, led_grid, spiral_order
from .model import FPM, ReconstructionResult, run_experiment
from .optimizers import simulated_annealing
from .solvers import fit_led_pose
from .visualize import save_summary, summary_figure

__all__ = [
    "DATA_ROOT",
    "FPM",
    "FT",
    "IntensityData",
    "OUTPUT_ROOT",
    "ReconstructionResult",
    "dataset_name",
    "fit_led_pose",
    "iFT",
    "illumination_wavevectors",
    "led_grid",
    "load_config",
    "load_intensity",
    "matlab_round",
    "output_dir",
    "run_experiment",
    "save_summary",
    "simulated_annealing",
    "spiral_order",
    "summary_figure",
]
