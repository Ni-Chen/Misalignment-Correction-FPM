"""The Fourier ptychographic microscope as a single configured object.

Everything that depends on a particular experiment -- the optics, the LED
array, the measurements and the solvers -- lives on :class:`FPM`. Stateless
mathematics (Fourier transforms, LED geometry, simulated annealing) lives in
the sibling modules and is imported here.

    from fpm import FPM, output_dir

    model = FPM.from_directory("data/usaf/edge_upper_right", method="mcfpm")
    result = model.run()
    result.save(output_dir("data/usaf/edge_upper_right", "mcfpm"), model=model, summary=True)

Conventions shared by every solver:

* A spectrum centre is the continuous pixel coordinate, 0-based, at which an
  LED's wavevector lands in the ``size_hr`` object spectrum:
  ``(size_hr - 1) / 2 + k / dk``. Its window is the ``size_lr`` block starting
  at ``round_half_away(centre - (size_lr - 1) / 2)``.
* The low-resolution field of one LED is
  ``iFT(scale * window * transfer_function * pupil * defocus_kernel)`` with
  ``scale = (size_lr / size_hr)**2``, and its modelled amplitude is ``scale``
  times the modulus of that field, so the modulus is replaced by
  ``sqrt(I) / scale``.
* The object spectrum starts as the zero-padded spectrum of the centre frame's
  *intensity*, and the pupil as ones.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch

from .config import DATA_ROOT, load_config
from .data import IntensityData, load_intensity
from .fourier import FT, iFT, round_half_away
from .illumination import illumination_wavevectors, led_grid, spiral_order
from .optimizers import simulated_annealing

# Every solver here works on one frame-sized array at a time, and MKL's FFT
# does not parallelise at that size: a 128 x 128 complex128 transform took
# 0.05 ms on 1-8 threads but 9 ms on all 24 logical cores of a laptop i9,
# which made a sweep 40x slower. Solvers on the CPU therefore run on at most
# this many threads; the setting is restored afterwards.
CPU_THREADS = 4


@contextmanager
def _limited_cpu_threads(device: torch.device):
    if device.type != "cpu" or torch.get_num_threads() <= CPU_THREADS:
        yield
        return
    previous = torch.get_num_threads()
    torch.set_num_threads(CPU_THREADS)
    try:
        yield
    finally:
        torch.set_num_threads(previous)


class FPM:
    """Config-driven Fourier ptychographic microscope.

    Measurements are held as intensity, exactly as stored on disk. The
    solvers take ``sqrt(intensity)`` themselves where the update needs an
    amplitude.
    """

    # Every solver's defaults; data/_template lists them per method.
    DEFAULTS = {
        "upscale": 4,         # object side length / measured side length
        "alpha": 1.0,         # object-update regulariser, delta_1 in the paper
        "beta": 1.0,          # pupil-update regulariser, delta_2 in the paper, which states 1000
        "seed": 0,
        "print_every": 10,    # progress of the long searches; 0 disables
    }

    def __init__(
        self,
        config: dict[str, Any],
        intensity_data: IntensityData,
        *,
        device: str | torch.device | None = None,
        dtype: torch.dtype = torch.float64,
    ) -> None:
        self.config = config
        self.method = str(config["method"])
        self.dataset_dir = Path(config["dataset_dir"])
        self.dtype = dtype
        self.complex_dtype = torch.complex128 if dtype == torch.float64 else torch.complex64
        self.device = torch.device(
            device if device is not None else ("cuda" if torch.cuda.is_available() else "cpu")
        )

        self.intensity_data = intensity_data
        self.intensity = torch.as_tensor(intensity_data.values, device=self.device, dtype=self.dtype)
        self.amplitude = torch.sqrt(self.intensity)

        optics = config["optics"]
        led = config["led"]
        model = config.get("model") or {}
        self.recon = {**self.DEFAULTS, **(config.get("reconstruction") or {})}

        # Optical system
        self.wavelength = float(optics["wavelength_m"])
        self.na = float(optics["na"])
        self.magnification = float(optics["magnification"])
        self.pixel_size_sensor = float(optics["pixel_size_sensor_m"])
        self.defocus = float(model.get("defocus_m", 0.0))
        self.upscale = int(self.recon["upscale"])

        # LED array, in millimetres
        self.array_size = int(led["array_size"])
        self.led_pitch_mm = float(led["pitch_m"]) * 1e3
        self.led_height_mm = float(led["height_m"]) * 1e3
        self.sample_offset_mm = (
            float(model.get("sample_offset_x_mm", 0.0)),
            float(model.get("sample_offset_y_mm", 0.0)),
        )
        self.led_grid = led_grid(self.array_size, self.led_pitch_mm)
        self.scan_order = spiral_order(self.array_size)
        if self.led_grid.shape[0] != self.n_frames:
            raise ValueError(
                f"LED count ({self.led_grid.shape[0]}) does not match "
                f"measurement count ({self.n_frames})."
            )

        self.transfer_function = self._transfer_function()
        self.defocus_kernel = self._defocus_kernel()
        self.last_result: ReconstructionResult | None = None

    @classmethod
    def from_directory(
        cls,
        dataset_dir: str | Path,
        method: str,
        *,
        device: str | torch.device | None = None,
        dtype: torch.dtype = torch.float64,
    ) -> "FPM":
        """Build a model from ``config.yaml`` plus the selected method config."""
        config = load_config(dataset_dir, method)
        return cls(config, load_intensity(config), device=device, dtype=dtype)

    def __repr__(self) -> str:
        return (
            f"FPM(dataset={self.dataset_dir.name!r}, method={self.method!r}, "
            f"frames={self.n_frames}, size={self.size_lr}->{self.size_hr}, "
            f"device={self.device})"
        )

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------

    @property
    def n_frames(self) -> int:
        """Number of recorded low-resolution frames."""
        return int(self.intensity.shape[0])

    @property
    def size_lr(self) -> int:
        """Side length of one recorded frame, in pixels."""
        return int(self.intensity.shape[-1])

    @property
    def size_hr(self) -> int:
        """Side length of the reconstructed object, in pixels."""
        return self.upscale * self.size_lr

    @property
    def pixel_size_object(self) -> float:
        """Sensor pixel pitch referred to the sample plane."""
        return self.pixel_size_sensor / self.magnification

    @property
    def k0(self) -> float:
        """Free-space wavenumber."""
        return 2 * np.pi / self.wavelength

    @property
    def dk(self) -> float:
        """Frequency step of the object spectrum, the same along both axes."""
        return 2 * np.pi / (self.pixel_size_object / self.upscale * self.size_hr)

    @property
    def centre_frame(self) -> int:
        """Raster index of the LED at the centre of the array, the first one the spiral visits."""
        return int(self.scan_order[0])

    @property
    def nominal_pose(self) -> np.ndarray:
        """``(x_mm, y_mm, theta_rad, height_mm)`` of the array as built: no shift, no rotation."""
        return np.array([0.0, 0.0, 0.0, self.led_height_mm])

    @property
    def _scale(self) -> float:
        """``(size_lr / size_hr)**2``, the energy scale between the object and one frame."""
        return (self.size_lr / self.size_hr) ** 2

    def _frequency_grid(self):
        k_max = np.pi / self.pixel_size_object
        k = np.linspace(-k_max, k_max, self.size_lr)
        return np.meshgrid(k, k)

    def _transfer_function(self) -> torch.Tensor:
        """Binary coherent transfer function of the objective (``CTF``)."""
        kx, ky = self._frequency_grid()
        support = (kx ** 2 + ky ** 2) < (self.na * self.k0) ** 2
        return torch.as_tensor(support, device=self.device, dtype=self.complex_dtype)

    def _defocus_kernel(self) -> torch.Tensor:
        """Angular-spectrum propagator over ``defocus_m``; ones in focus."""
        kx, ky = self._frequency_grid()
        kz = np.sqrt((self.k0 ** 2 - kx ** 2 - ky ** 2).astype(np.complex128))
        kernel = np.exp(1j * self.defocus * kz.real) * np.exp(-self.defocus * kz.imag)
        return torch.as_tensor(kernel, device=self.device, dtype=self.complex_dtype)

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def spectrum_centres(self, pose=None) -> np.ndarray:
        """``(N, 2)`` centres ``(x, y)`` of every LED's window for an array pose.

        ``pose`` defaults to :attr:`nominal_pose`. Each centre is the LED's
        wavevector divided by :attr:`dk`, counted from the spectrum origin
        ``(size_hr - 1) / 2``.
        """
        pose = self.nominal_pose if pose is None else np.asarray(pose, dtype=np.float64)
        kx, ky = illumination_wavevectors(self.led_grid, pose, self.sample_offset_mm, self.k0)
        origin = (self.size_hr - 1) / 2
        return np.stack([origin + kx / self.dk, origin + ky / self.dk], axis=1)

    def window(self, centre) -> tuple[int, int]:
        """Top-left corner ``(y0, x0)`` of the window centred on ``centre = (x, y)``."""
        half = (self.size_lr - 1) / 2
        x0, y0 = (int(value) for value in round_half_away(np.asarray(centre) - half))
        if not (0 <= x0 <= self.size_hr - self.size_lr and 0 <= y0 <= self.size_hr - self.size_lr):
            raise ValueError(
                f"Spectrum window at centre {tuple(np.round(centre, 2))} falls outside the "
                f"{self.size_hr} px object spectrum; the LED pose or sample offset is out of range."
            )
        return y0, x0

    @staticmethod
    def pose_bounds(settings: dict[str, Any], nominal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Absolute ``(lower, upper)`` pose bounds from bounds stated relative to ``nominal``.

        ``settings`` maps ``x_mm``, ``y_mm``, ``rotation_deg`` and ``height_mm``
        to ``[low, high]`` offsets; the rotation is converted to radians.
        """
        keys = ("x_mm", "y_mm", "rotation_deg", "height_mm")
        missing = [key for key in keys if key not in settings]
        if missing:
            raise KeyError(f"pose_bounds is missing {missing}")
        offsets = np.array([settings[key] for key in keys], dtype=np.float64)
        offsets[2] = np.deg2rad(offsets[2])
        return nominal + offsets[:, 0], nominal + offsets[:, 1]

    # ------------------------------------------------------------------
    # The ptychographic engine
    # ------------------------------------------------------------------

    def initial_spectrum(self) -> torch.Tensor:
        """Zero-padded spectrum of the centre frame's intensity."""
        pad = (self.size_hr - self.size_lr) // 2
        spectrum = FT(self.intensity[self.centre_frame].to(self.complex_dtype))
        return torch.nn.functional.pad(spectrum, (pad, pad, pad, pad))

    def initial_pupil(self) -> torch.Tensor:
        return torch.ones((self.size_lr, self.size_lr), device=self.device, dtype=self.complex_dtype)

    def _update(self, spectrum, pupil, centre, frame: int, cost: str):
        """One LED's object and pupil update; ``spectrum`` is changed in place.

        Returns the new pupil and the error of the prediction made before the
        update, summed over the frame: ``(scale |psi| - I)**2`` for
        ``cost="amplitude_vs_intensity"`` or ``((scale |psi|)**2 - I)**2`` for
        ``"intensity"`` (Eq. 15 of the paper).
        """
        alpha = float(self.recon["alpha"])
        beta = float(self.recon["beta"])
        y0, x0 = self.window(centre)
        size = self.size_lr
        patch = spectrum[y0:y0 + size, x0:x0 + size]

        predicted_spectrum = self._scale * patch * self.transfer_function * pupil * self.defocus_kernel
        predicted_field = iFT(predicted_spectrum)
        predicted_amplitude = self._scale * predicted_field.abs()
        if cost == "amplitude_vs_intensity":
            error = (predicted_amplitude - self.intensity[frame]).square().sum()
        elif cost == "intensity":
            error = (predicted_amplitude.square() - self.intensity[frame]).square().sum()
        else:
            raise ValueError(f"cost must be 'amplitude_vs_intensity' or 'intensity', got {cost!r}")

        # Replace the modulus by the measurement, keep the phase.
        corrected_field = self.upscale ** 2 * self.amplitude[frame] * torch.exp(1j * predicted_field.angle())
        difference = FT(corrected_field) - predicted_spectrum

        pupil_magnitude = pupil.abs()
        patch = patch + (
            1 / pupil_magnitude.max() * pupil_magnitude * pupil.conj() * difference
            / (pupil_magnitude.square() + alpha)
        )
        spectrum[y0:y0 + size, x0:x0 + size] = patch

        patch_magnitude = patch.abs()
        pupil = pupil + (
            1 / patch_magnitude.max() * patch_magnitude * patch.conj() * difference
            / (patch_magnitude.square() + beta) * self.transfer_function
        )
        return pupil, error

    def position_cost(self, spectrum, pupil, centre, frame: int) -> float:
        """Intensity misfit of one frame for a trial window centre."""
        y0, x0 = self.window(centre)
        size = self.size_lr
        patch = spectrum[y0:y0 + size, x0:x0 + size]
        field = iFT(self._scale * self.transfer_function * self.defocus_kernel * pupil * patch)
        predicted = (self._scale * field).abs().square()
        return float((predicted - self.intensity[frame]).square().sum())

    def search_position(self, spectrum, pupil, centre, frame: int, radius: float, seed) -> np.ndarray:
        """Anneal one LED's window centre within ``radius`` pixels of ``centre``."""
        centre = np.asarray(centre, dtype=np.float64)
        result = simulated_annealing(
            lambda trial: self.position_cost(spectrum, pupil, trial, frame),
            centre,
            centre - radius,
            centre + radius,
            function_tolerance=float(self.recon.get("search_function_tolerance", 1e-3)),
            stall_iterations=self.recon.get("search_stall_iterations"),
            max_function_evaluations=self.recon.get("search_max_function_evaluations"),
            seed=seed,
        )
        return result.x

    def sweep(self, spectrum, pupil, centres, order, *, cost: str = "intensity",
              search_radius: float | None = None, rng: np.random.Generator | None = None):
        """Update every LED of ``order`` once.

        With ``search_radius`` each LED's centre is first annealed and written
        back into ``centres``. Returns the pupil and the summed error.
        """
        error = torch.zeros((), device=self.device, dtype=self.dtype)
        for frame in order:
            if search_radius is not None:
                seed = None if rng is None else int(rng.integers(2 ** 32))
                centres[frame] = self.search_position(spectrum, pupil, centres[frame], frame,
                                                      search_radius, seed)
            pupil, frame_error = self._update(spectrum, pupil, centres[frame], frame, cost)
            error = error + frame_error
        return pupil, float(error)

    def reconstruct(self, centres, iterations: int, *, order=None, cost: str = "intensity",
                    plain_sweeps: int | None = None, search_radius: float | None = None,
                    rng: np.random.Generator | None = None):
        """Sequential reconstruction from a fresh start, the loop every solver is built on.

        The first ``plain_sweeps`` sweeps keep ``centres`` fixed; every later
        one anneals each LED's centre within ``search_radius`` first.
        Returns the object spectrum, the pupil, the centres actually used and
        the error of every sweep.
        """
        order = self.scan_order if order is None else order
        centres = np.array(centres, dtype=np.float64)
        errors = []
        with _limited_cpu_threads(self.device):
            spectrum = self.initial_spectrum()
            pupil = self.initial_pupil()
            for iteration in range(iterations):
                searching = plain_sweeps is not None and iteration >= plain_sweeps
                pupil, error = self.sweep(spectrum, pupil, centres, order, cost=cost,
                                          search_radius=search_radius if searching else None, rng=rng)
                errors.append(error)
        return spectrum, pupil, centres, errors

    # ------------------------------------------------------------------
    # Solvers
    # ------------------------------------------------------------------

    def run(self) -> ReconstructionResult:
        """Run the solver named by ``config['method']``."""
        with _limited_cpu_threads(self.device):
            return getattr(self, self.method)()

    def pie(self) -> ReconstructionResult:
        """Conventional FPM: the LED array at its nominal pose, no correction.

        Each centre is rounded to a whole pixel before its window is cut,
        unlike in the misalignment solvers. With an even object size the
        nominal centres sit on half pixels, so the rounding moves about half
        of the windows by a pixel (``test_rounding_the_centres_costs_pie_accuracy``).

        ``result.loss`` holds the intensity misfit of every sweep (Eq. 15 of
        the paper), measured before each LED's update.
        """
        iterations = int(self.recon.get("iterations", 9))
        centres = round_half_away(self.spectrum_centres())
        spectrum, pupil, centres, errors = self.reconstruct(centres, iterations)
        return self._result("pie", spectrum, pupil, centres, errors, pose=self.nominal_pose)

    def sa(self) -> ReconstructionResult:
        """Conventional simulated-annealing position correction.

        After ``plain_sweeps`` ordinary sweeps from the rounded nominal
        centres, every sweep anneals each LED's centre independently within
        ``search_radius_px`` before updating with it. The LEDs are not tied
        to an array model, which is what lets them drift apart.
        """
        iterations = int(self.recon.get("iterations", 10))
        plain_sweeps = int(self.recon.get("plain_sweeps", 9))
        radius = float(self.recon.get("search_radius_px", 10.0))
        rng = np.random.default_rng(int(self.recon["seed"]))
        centres = round_half_away(self.spectrum_centres())
        spectrum, pupil, centres, errors = self.reconstruct(
            centres, iterations, plain_sweeps=plain_sweeps, search_radius=radius, rng=rng,
        )
        return self._result("sa", spectrum, pupil, centres, errors, pose=None)

    def pcfpm(self) -> ReconstructionResult:
        """SA plus a global array model. See :func:`fpm.solvers.pcfpm`."""
        from .solvers import pcfpm

        return pcfpm(self)

    def mcfpm(self) -> ReconstructionResult:
        """Global-shift misalignment correction. See :func:`fpm.solvers.mcfpm`."""
        from .solvers import mcfpm

        return mcfpm(self)

    def _result(self, method: str, spectrum, pupil, centres, errors, *, pose,
                **extra: Any) -> ReconstructionResult:
        result = ReconstructionResult(
            method=method,
            object=iFT(spectrum).cpu(),
            pupil=pupil.cpu(),
            spectrum_centres=np.asarray(centres),
            initial_spectrum_centres=self.spectrum_centres(),
            led_pose=None if pose is None else _pose_for_humans(pose),
            loss=list(errors),
            metadata=self._metadata(),
            **extra,
        )
        self.last_result = result
        return result

    def _metadata(self) -> dict[str, Any]:
        paths = self.config.get("config_paths") or {}
        return {
            "dataset": str(self.dataset_dir),
            "data_file": str(self.intensity_data.source_path),
            "config_common": str(paths.get("common")),
            "config_method": str(paths.get("method")),
            "reconstruction": self.recon,
            "sample_offset_mm": list(self.sample_offset_mm),
            "upscale": self.upscale,
            "dtype": str(self.dtype),
            "device": str(self.device),
            "pose_convention": "x_mm, y_mm, rotation_deg, height_mm",
        }

    def summary_figure(self, result: ReconstructionResult | None = None, **kwargs: Any):
        """Overview figure of ``result`` (default: the last one this model produced)."""
        from .visualize import summary_figure

        result = result if result is not None else self.last_result
        if result is None:
            raise RuntimeError("No reconstruction to show; call run() first.")
        return summary_figure(result, self, **kwargs)


def _pose_for_humans(pose) -> np.ndarray:
    """The pose with its rotation in degrees, as results report it."""
    pose = np.array(pose, dtype=np.float64)
    pose[2] = np.rad2deg(pose[2])
    return pose


@dataclass
class ReconstructionResult:
    """Object, pupil, LED geometry and traces returned by a solver."""

    method: str
    object: torch.Tensor                        # (H, W) complex high-resolution field
    pupil: torch.Tensor                         # (h, w) complex pupil
    spectrum_centres: np.ndarray                # (N, 2) final window centres (x, y), 0-based pixels
    initial_spectrum_centres: np.ndarray        # (N, 2) centres of the nominal array
    led_pose: np.ndarray | None = None          # (x_mm, y_mm, rotation_deg, height_mm); None for SA
    led_pose_history: list[list[float]] = field(default_factory=list)
    loss: list[float] = field(default_factory=list)             # per sweep of the final reconstruction
    anneal_history: list[float] = field(default_factory=list)   # mcFPM: best pose cost per annealing iteration
    metadata: dict[str, Any] = field(default_factory=dict)

    def figure(self, model: "FPM | None" = None, **kwargs: Any):
        """Build the multi-panel summary figure for this reconstruction."""
        from .visualize import summary_figure

        return summary_figure(self, model, **kwargs)

    def save(self, output_dir: str | Path, *, model: "FPM | None" = None,
             summary: bool = False) -> Path:
        """Write ``reconstruction.npz`` and ``metadata.json`` into ``output_dir``.

        With ``summary=True`` a ``summary.png`` overview is rendered as well.
        Refuses to write under ``data/``, which holds acquisitions only; use
        ``fpm.config.output_dir(dataset_dir, method)``.
        """
        output_path = Path(output_dir)
        resolved = output_path.expanduser().resolve()
        try:
            resolved.relative_to(DATA_ROOT.resolve())
        except ValueError:
            pass
        else:
            raise ValueError(
                f"Refusing to save a reconstruction under {DATA_ROOT} ({resolved}): "
                "that tree holds acquisitions only. Use "
                "fpm.config.output_dir(dataset_dir, method) instead."
            )
        output_path.mkdir(parents=True, exist_ok=True)
        arrays: dict[str, np.ndarray] = {
            "object": self.object.numpy(),
            "pupil": self.pupil.numpy(),
            "spectrum_centres": self.spectrum_centres,
            "initial_spectrum_centres": self.initial_spectrum_centres,
            "loss": np.asarray(self.loss, dtype=np.float64),
            "anneal_history": np.asarray(self.anneal_history, dtype=np.float64),
        }
        if self.led_pose is not None:
            arrays["led_pose"] = self.led_pose
        if self.led_pose_history:
            arrays["led_pose_history"] = np.asarray(self.led_pose_history, dtype=np.float64)
        np.savez_compressed(output_path / "reconstruction.npz", **arrays)
        metadata = {"method": self.method, **self.metadata}
        if self.led_pose is not None:
            metadata["led_pose"] = self.led_pose.tolist()
        (output_path / "metadata.json").write_text(
            json.dumps(metadata, indent=2, default=str), encoding="utf-8",
        )
        if summary:
            from .visualize import save_summary

            save_summary(self, output_path, model)
        return output_path


def run_experiment(
    dataset_dir: str | Path,
    method: str,
    *,
    device: str | torch.device | None = None,
    output_dir: str | Path | None = None,
    summary: bool = True,
) -> ReconstructionResult:
    """Load a dataset directory, run its configured solver, optionally save."""
    model = FPM.from_directory(dataset_dir, method, device=device)
    result = model.run()
    if output_dir is not None:
        result.save(output_dir, model=model, summary=summary)
    return result
