"""Compare saved reconstructions of every dataset side by side.

    python tools/summarize_results.py                        # pie against mcfpm
    python tools/summarize_results.py --methods pie pcfpm
    python tools/summarize_results.py --wide --refit 9 --output docs/img/comparison_pie_mcfpm.png

or from Python as ``summarize_results.main(["--methods", "pie", "mcfpm"])``,
which also returns the figure.

Reads ``reconstruction.npz`` and ``metadata.json`` from
``notebooks/output/<dataset>/<method>/`` -- what is on disk, not what a
config says should happen -- and writes one figure with the amplitude and
phase of each method for every dataset, a row per dataset or, with
``--wide``, a column per dataset, to ``notebooks/output/comparison_<methods>.png``
unless ``--output`` says otherwise. A table of the patch offsets, recovered
LED poses and final losses is printed. The losses are each solver's
last-sweep intensity misfit and are not comparable between methods that run
different numbers of sweeps.

``--refit N`` makes the comparison that is: for every saved result with an
annealed LED pose it reconstructs N sweeps with the same engine, once with
the uncorrected array and once with the annealed pose, and reports both
misfits. It also converts the annealed shift into the patch position it is
equivalent to -- a shift of the array and of the patch enter the wavevector
only through their difference -- and divides that by the position the patch
was cut at (``config.yaml``).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fpm import FPM
from fpm.config import DATA_ROOT, METHOD_FILES, OUTPUT_ROOT, dataset_name, iter_dataset_dirs, output_dir

FONT_SIZE = 8


def load(dataset: Path, method: str):
    directory = output_dir(dataset, method)
    if not (directory / "reconstruction.npz").is_file():
        return None
    with np.load(directory / "reconstruction.npz") as arrays:
        result = {key: arrays[key] for key in arrays.files}
    result["metadata"] = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    return result


def without_global_phase(field: np.ndarray) -> np.ndarray:
    """Rotate the field so its mean is real; the absolute phase is not measurable."""
    return field * np.exp(-1j * np.angle(field.mean()))


def draw(ax, result, part: str) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    if result is None:
        ax.text(0.5, 0.5, "not run", ha="center", va="center", fontsize=FONT_SIZE, color="0.5")
    elif part == "amplitude":
        ax.imshow(np.abs(result["object"]), cmap="gray")
    else:
        phase = np.angle(without_global_phase(result["object"]))
        limit = max(np.percentile(np.abs(phase), 99.5), 0.05)
        ax.imshow(phase, cmap="viridis", vmin=-limit, vmax=limit)


def figure_of(datasets, results, methods, wide: bool):
    panels = [(method, part) for method in methods for part in ("amplitude", "phase")]
    if wide:
        panels = [(method, part) for part in ("amplitude", "phase") for method in methods]
        figure, axes = plt.subplots(len(panels), len(datasets),
                                    figsize=(1.75 * len(datasets), 1.75 * len(panels)), squeeze=False)
        for column, dataset in enumerate(datasets):
            for row, (method, part) in enumerate(panels):
                ax = axes[row, column]
                draw(ax, results[dataset][method], part)
                if row == 0:
                    ax.set_title(dataset_name(dataset).replace("/", "/\n"), fontsize=FONT_SIZE)
                if column == 0:
                    ax.set_ylabel(f"{method}\n{part}", fontsize=FONT_SIZE + 1)
    else:
        figure, axes = plt.subplots(len(datasets), len(panels),
                                    figsize=(2.3 * len(panels), 2.3 * len(datasets)), squeeze=False)
        for row, dataset in enumerate(datasets):
            for column, (method, part) in enumerate(panels):
                ax = axes[row, column]
                draw(ax, results[dataset][method], part)
                if row == 0:
                    ax.set_title(f"{method} {part}", fontsize=FONT_SIZE + 1)
                if column == 0:
                    ax.set_ylabel(dataset_name(dataset), fontsize=FONT_SIZE)
    figure.tight_layout()
    return figure


def refit(datasets, results, methods, sweeps: int) -> None:
    """Misfit of the uncorrected and the annealed geometry after the same reconstruction."""
    print(f"\nSame engine, {sweeps} sweeps, uncorrected array against the annealed pose:\n")
    print("| dataset | method | patch cut at (x, y) mm | annealed shift (x, y) mm | "
          "misfit uncorrected | misfit annealed | ratio | equivalent patch position / cut position (x, y) |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for dataset in datasets:
        model_config = yaml.safe_load((dataset / "config.yaml").read_text(encoding="utf-8")).get("model") or {}
        cut = np.array([model_config.get("sample_offset_x_mm", 0.0), model_config.get("sample_offset_y_mm", 0.0)])
        for method in methods:
            result = results[dataset][method]
            if result is None or "led_pose" not in result or not np.any(result["led_pose"][:2]):
                continue
            pose = np.array(result["led_pose"], dtype=np.float64)
            pose[2] = np.deg2rad(pose[2])
            model = FPM.from_directory(dataset, method)
            uncorrected = model.reconstruct(model.spectrum_centres(), sweeps)[3][-1]
            annealed = model.reconstruct(model.spectrum_centres(pose), sweeps)[3][-1]
            equivalent = np.array(model.sample_offset_mm) - pose[:2]
            ratio = "" if not np.all(np.abs(cut) > 1e-9) else ", ".join(f"{v:.2f}" for v in equivalent / cut)
            print(f"| {dataset_name(dataset)} | {method} | {cut[0]:+.3f}, {cut[1]:+.3f} | "
                  f"{pose[0]:+.3f}, {pose[1]:+.3f} | {uncorrected:.3g} | {annealed:.3g} | "
                  f"{uncorrected / annealed:.2f} | {ratio} |", flush=True)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--methods", nargs="+", choices=tuple(METHOD_FILES), default=["pie", "mcfpm"])
    parser.add_argument("--wide", action="store_true", help="one column per dataset instead of one row")
    parser.add_argument("--output", type=Path, default=None, help="figure path (default notebooks/output/)")
    parser.add_argument("--refit", type=int, default=0, metavar="N",
                        help="also compare N-sweep reconstructions without and with the annealed pose")
    args = parser.parse_args(argv)

    datasets = list(iter_dataset_dirs(DATA_ROOT))
    results = {dataset: {method: load(dataset, method) for method in args.methods} for dataset in datasets}

    print("| dataset | method | patch offset (x, y) mm | LED pose (x mm, y mm, deg, h mm) | final loss |")
    print("| --- | --- | --- | --- | --- |")
    for dataset in datasets:
        for method, result in results[dataset].items():
            if result is None:
                continue
            offset = ", ".join(f"{value:.4g}" for value in result["metadata"].get("sample_offset_mm"))
            pose = result.get("led_pose")
            pose = "" if pose is None else ", ".join(f"{value:.3f}" for value in pose)
            loss = result["loss"][-1] if result["loss"].size else float("nan")
            print(f"| {dataset_name(dataset)} | {method} | {offset} | {pose} | {loss:.3g} |")

    figure = figure_of(datasets, results, args.methods, args.wide)
    target = args.output or OUTPUT_ROOT / f"comparison_{'_'.join(args.methods)}.png"
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=110)
    print(f"\nWrote {target}")

    if args.refit:
        refit(datasets, results, args.methods, args.refit)
    return figure


if __name__ == "__main__":
    matplotlib.use("Agg")  # file output only; a notebook keeps its own backend
    main()
