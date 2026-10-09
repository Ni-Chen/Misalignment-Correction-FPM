"""Command-line entry point: ``python -m fpm <dataset_dir> --method mcfpm``."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import METHOD_FILES
from .config import output_dir as default_output_dir
from .model import run_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a config-driven FPM reconstruction.")
    parser.add_argument("dataset_dir", type=Path, help="Directory containing data and config files")
    parser.add_argument("--method", choices=tuple(METHOD_FILES), required=True)
    parser.add_argument("--device", default=None, help="torch device, for example cuda:0 or cpu")
    parser.add_argument("--output", type=Path, default=None,
                        help="override the default notebooks/output/<dataset>/<method>")
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="skip the summary.png overview figure",
    )
    args = parser.parse_args()

    output = args.output or default_output_dir(args.dataset_dir, args.method)
    result = run_experiment(
        args.dataset_dir,
        args.method,
        device=args.device,
        output_dir=output,
        summary=not args.no_summary,
    )
    print(f"Saved {result.method} reconstruction to {output.resolve()}")


if __name__ == "__main__":
    main()
