"""Result figures.

``summary_figure`` puts one reconstruction on a page: the recorded centre
frame, amplitude and phase, the object spectrum, the pupil, where the LED
windows ended up, and the convergence traces. The panels are separate
functions so they can be reused for paper figures.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

AMPLITUDE_CMAP = "gray"
PHASE_CMAP = "viridis"
SPECTRUM_CMAP = "magma"
FONT_SIZE = 9


def _to_numpy(value) -> np.ndarray:
    if hasattr(value, "detach"):
        value = value.detach().cpu().numpy()
    return np.asarray(value)


def _image(ax, values, title: str, *, cmap: str = AMPLITUDE_CMAP) -> None:
    image = ax.imshow(values, cmap=cmap)
    ax.set_title(title, fontsize=FONT_SIZE)
    ax.set_xticks([])
    ax.set_yticks([])
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.ax.tick_params(labelsize=FONT_SIZE - 1)


def _blank(ax, message: str) -> None:
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=FONT_SIZE, color="0.5")
    ax.set_axis_off()


def plot_amplitude(ax, object_field) -> None:
    _image(ax, np.abs(_to_numpy(object_field)), "Amplitude")


def plot_phase(ax, object_field) -> None:
    _image(ax, np.angle(_to_numpy(object_field)), "Phase (rad)", cmap=PHASE_CMAP)


def plot_spectrum(ax, object_field) -> None:
    spectrum = np.fft.fftshift(np.fft.fft2(_to_numpy(object_field)))
    _image(ax, np.log10(np.abs(spectrum) + 1e-12), "log10 |spectrum|", cmap=SPECTRUM_CMAP)


def plot_pupil(ax_amplitude, ax_phase, pupil) -> None:
    pupil = _to_numpy(pupil)
    _image(ax_amplitude, np.abs(pupil), "Pupil amplitude")
    _image(ax_phase, np.angle(pupil), "Pupil phase (rad)", cmap=PHASE_CMAP)


def plot_spectrum_centres(ax, result) -> None:
    """Nominal against final window centres."""
    nominal = np.asarray(result.initial_spectrum_centres)
    final = np.asarray(result.spectrum_centres)
    ax.scatter(nominal[:, 0], nominal[:, 1], s=12, marker="x", color="k", label="nominal")
    ax.scatter(final[:, 0], final[:, 1], s=14, facecolors="none", edgecolors="tab:red", label="corrected")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    title = "Spectrum centres (px)"
    if result.led_pose is not None:
        x, y, rotation, height = result.led_pose
        title += f"\nshift ({x:.3f}, {y:.3f}) mm, {rotation:.3f} deg, h {height:.2f} mm"
    ax.set_title(title, fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE - 1)
    ax.legend(fontsize=FONT_SIZE - 1, loc="upper right")


def plot_trace(ax, values, title: str, xlabel: str, *, first: int = 0) -> None:
    """``values`` against their index, counted from ``first``."""
    if len(values) == 0:
        _blank(ax, f"no {title.lower()}")
        return
    ax.semilogy(np.arange(first, first + len(values)), values, marker=".")
    ax.set_title(title, fontsize=FONT_SIZE)
    ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE - 1)


def summary_figure(result, model=None, *, title: str | None = None):
    """One page: centre frame, amplitude, phase, spectrum, pupil, centres and traces."""
    figure, axes = plt.subplots(2, 4, figsize=(17, 8.5))
    axes = axes.ravel()
    if model is not None:
        _image(axes[0], _to_numpy(model.intensity[model.centre_frame]), "Centre frame (intensity)")
    else:
        _blank(axes[0], "no measurements")
    plot_amplitude(axes[1], result.object)
    plot_phase(axes[2], result.object)
    plot_spectrum(axes[3], result.object)
    plot_pupil(axes[4], axes[5], result.pupil)
    plot_spectrum_centres(axes[6], result)
    if result.anneal_history:
        # Entry 0 is the cost of the starting pose, before any annealing step.
        plot_trace(axes[7], result.anneal_history, "Annealing cost (best)", "annealing iteration")
    else:
        # One sweep updates every LED once; sweeps are numbered from 1.
        plot_trace(axes[7], result.loss, "Intensity misfit", "sweep over all LEDs", first=1)
    figure.suptitle(title or result.method, fontsize=FONT_SIZE + 2)
    figure.tight_layout()
    return figure


def save_summary(result, output_dir: str | Path, model=None, *, dpi: int = 150,
                 filename: str = "summary.png", **kwargs: Any) -> Path:
    """Render :func:`summary_figure` into ``output_dir/filename``."""
    target = Path(output_dir) / filename
    figure = summary_figure(result, model, **kwargs)
    figure.savefig(target, dpi=dpi)
    plt.close(figure)
    return target
