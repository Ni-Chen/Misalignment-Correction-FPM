"""Configuration loading, the datasets in ``data/``, and saving results."""

from __future__ import annotations

import json

import numpy as np
import pytest

from fpm import FPM, load_config, output_dir
from fpm.config import DATA_ROOT, METHOD_FILES, iter_dataset_dirs
from fpm.data import load_intensity
from support import ROOT, cpu_model, simulated_dataset, write_dataset

DATASETS = list(iter_dataset_dirs(DATA_ROOT))


def test_missing_method_file_takes_the_template_defaults(tmp_path):
    dataset = write_dataset(tmp_path, np.ones((9, 8, 8)), led={"array_size": 3})
    config = load_config(dataset, "mcfpm")
    assert config["reconstruction"] == {}
    assert config["config_paths"]["method"] is None
    model = FPM(config, load_intensity(config), device="cpu")
    assert model.recon["upscale"] == 4 and model.recon["beta"] == 1.0


def test_method_file_can_override_the_patch_position(tmp_path):
    dataset = write_dataset(tmp_path, np.ones((9, 8, 8)), led={"array_size": 3},
                            model={"sample_offset_x_mm": 0.1, "sample_offset_y_mm": 0.2},
                            pie={"iterations": 3, "model": {"sample_offset_x_mm": 2.0}})
    config = load_config(dataset, "pie")
    assert config["model"] == {"sample_offset_x_mm": 2.0, "sample_offset_y_mm": 0.2}
    assert config["reconstruction"] == {"iterations": 3}
    assert load_config(dataset, "sa")["model"]["sample_offset_x_mm"] == 0.1


def test_unknown_method_is_rejected(tmp_path):
    dataset = write_dataset(tmp_path, np.ones((9, 8, 8)), led={"array_size": 3})
    with pytest.raises(ValueError, match="Unknown reconstruction method"):
        load_config(dataset, "ad")


def test_template_lists_every_method():
    template = DATA_ROOT / "_template"
    assert sorted(path.name for path in template.glob("config_*.yaml")) == sorted(METHOD_FILES.values())


def test_every_dataset_is_found():
    names = sorted(path.relative_to(DATA_ROOT).as_posix() for path in DATASETS)
    assert names == ["root/center", "root/edge_upper_left", "root/edge_upper_right",
                     "root/inner_upper_left", "root/inner_upper_left_far", "root/middle_upper_left",
                     "usaf/center", "usaf/edge_upper_right"]


@pytest.mark.parametrize("dataset", DATASETS, ids=lambda path: path.relative_to(DATA_ROOT).as_posix())
@pytest.mark.parametrize("method", sorted(METHOD_FILES))
def test_every_dataset_loads_and_its_windows_fit(dataset, method):
    config = load_config(dataset, method)
    if not (dataset / config["data"]["file"]).is_file():
        pytest.skip("measurement file not present")
    model = FPM(config, load_intensity(config), device="cpu")
    assert model.n_frames == 289
    assert model.intensity.min() >= 0
    for centre in model.spectrum_centres():
        model.window(centre)


def test_save_writes_arrays_metadata_and_summary(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path / "data", pose=None, pie={"iterations": 1}), "pie")
    result = model.run()
    target = result.save(tmp_path / "out", model=model, summary=True)
    with np.load(target / "reconstruction.npz") as arrays:
        assert arrays["object"].shape == (model.size_hr, model.size_hr)
        assert arrays["spectrum_centres"].shape == (289, 2)
        assert arrays["led_pose"].tolist() == [0.0, 0.0, 0.0, 113.5]
    metadata = json.loads((target / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["method"] == "pie" and metadata["reconstruction"]["iterations"] == 1
    assert (target / "summary.png").stat().st_size > 0


def test_results_are_never_written_under_data(tmp_path):
    model = cpu_model(simulated_dataset(tmp_path, pose=None, pie={"iterations": 1}), "pie")
    result = model.run()
    with pytest.raises(ValueError, match="acquisitions only"):
        result.save(DATA_ROOT / "usaf" / "edge_upper_right" / "result")
    assert output_dir(DATA_ROOT / "root" / "edge_upper_left", "pcfpm") == (
        ROOT / "notebooks" / "output" / "root" / "edge_upper_left" / "pcfpm"
    )
