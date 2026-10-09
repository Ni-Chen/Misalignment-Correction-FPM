"""Reconstruct every dataset with one or more methods and save results plus summary figures.

    python tools/run_all.py                                  # pie, then mcfpm, on every dataset
    python tools/run_all.py --methods pcfpm --only root      # one method on one group
    python tools/run_all.py --skip-existing                  # fill the gaps only

or from Python (a notebook) as ``run_all.main(["--methods", "pie"])``.

Methods run in the order given, each over all datasets before the next
starts. Each run is isolated and a failure is reported rather than fatal, so
one bad configuration does not cost the whole batch. Everything printed is
also written to ``<output_root>/run_all/<timestamp>.log`` and the outcome of
every run to ``<output_root>/run_all/<timestamp>.json`` beside it.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import sys
import time
import traceback

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fpm import FPM
from fpm.config import DATA_ROOT, METHOD_FILES, OUTPUT_ROOT, dataset_name, iter_dataset_dirs, output_dir


class Tee:
    """Write to the console and to a log file at the same time."""

    def __init__(self, stream, log) -> None:
        self._stream = stream
        self._log = log

    def write(self, text: str) -> int:
        self._stream.write(text)
        self._log.write(text)
        return len(text)

    def flush(self) -> None:
        self._stream.flush()
        self._log.flush()

    def __getattr__(self, name: str):
        return getattr(self._stream, name)


def selected_datasets(only: list[str] | None) -> list[Path]:
    """Every dataset, or those whose name equals or lies under one of ``only`` (``root`` is the whole group)."""
    datasets = list(iter_dataset_dirs(DATA_ROOT))
    if not only:
        return datasets
    keep = [path for path in datasets
            if any(dataset_name(path) == name or dataset_name(path).startswith(name.rstrip("/") + "/")
                   for name in only)]
    if not keep:
        raise SystemExit(f"No dataset matches {only}")
    return keep


def run_one(dataset: Path, method: str, device: str | None) -> dict:
    started = time.perf_counter()
    model = FPM.from_directory(dataset, method, device=device)
    result = model.run()
    target = result.save(output_dir(dataset, method), model=model, summary=True)
    record = {
        "seconds": round(time.perf_counter() - started, 1),
        "output": str(target),
        "final_loss": result.loss[-1] if result.loss else None,
        "sample_offset_mm": list(model.sample_offset_mm),
    }
    if result.led_pose is not None:
        record["led_pose"] = [round(float(value), 4) for value in result.led_pose]
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--methods", nargs="+", choices=tuple(METHOD_FILES), default=["pie", "mcfpm"])
    parser.add_argument("--only", nargs="+", default=None, help="dataset names or groups, e.g. root usaf/center")
    parser.add_argument("--skip-existing", action="store_true", help="skip runs whose reconstruction.npz exists")
    parser.add_argument("--device", default=None, help="torch device, for example cuda:0 or cpu")
    args = parser.parse_args(argv)

    log_dir = OUTPUT_ROOT / "run_all"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    log = (log_dir / f"{stamp}.log").open("w", encoding="utf-8")
    console = sys.stdout
    sys.stdout = Tee(console, log)
    try:
        records = run(args, log_dir, stamp)
    finally:
        sys.stdout = console
        log.close()
    return 1 if any(record["status"] == "failed" for record in records) else 0


def run(args, log_dir: Path, stamp: str) -> list[dict]:
    datasets = selected_datasets(args.only)
    print(f"{len(datasets)} datasets x {args.methods}", flush=True)
    records = []
    for method in args.methods:
        for dataset in datasets:
            name = dataset_name(dataset)
            record = {"dataset": name, "method": method}
            if args.skip_existing and (output_dir(dataset, method) / "reconstruction.npz").is_file():
                record["status"] = "skipped"
                print(f"[skip] {name} {method}", flush=True)
            else:
                print(f"[run ] {name} {method}", flush=True)
                try:
                    record.update(run_one(dataset, method, args.device), status="ok")
                    print(f"[done] {name} {method} in {record['seconds']} s", flush=True)
                except Exception as error:
                    record.update(status="failed", error=f"{type(error).__name__}: {error}")
                    traceback.print_exc(file=sys.stdout)
                    print(f"[fail] {name} {method}: {record['error']}", flush=True)
            records.append(record)
            (log_dir / f"{stamp}.json").write_text(json.dumps(records, indent=2), encoding="utf-8")

    print("\ndataset                         method  status   seconds  final loss  LED pose (x mm, y mm, deg, h mm)")
    for record in records:
        loss = record.get("final_loss")
        print(f"{record['dataset']:31s} {record['method']:7s} {record['status']:8s} "
              f"{record.get('seconds', ''):>7}  {'' if loss is None else f'{loss:.4g}':>10}  "
              f"{record.get('led_pose', '')}", flush=True)
    return records


if __name__ == "__main__":
    raise SystemExit(main())
