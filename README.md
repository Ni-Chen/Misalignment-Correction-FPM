
# Fast and robust misalignment correction of Fourier ptychographic microscopy 

[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

Implementation of "Fast and robust misalignment correction of Fourier ptychographic microscopy for full field of view reconstruction" by Ao Zhou, Wei Wang, [Ni Chen*](https://ni-chen.github.io/), Edmund Y. Lam, Byoungho Lee, and Guohai Situ [[Paper]](https://doi.org/10.1364/OE.26.023661) [[Project page]](https://ni-chen.github.io/Misalignment-Correction-FPM/)

The reconstruction is a Python package, `src/fpm`. The original MATLAB code is archived.

## Installation

```powershell
uv venv
uv pip install -e ".[test]"
```

or `pip install -e ".[test]"` in any environment with Python 3.10 or newer. PyTorch runs on the GPU when one is available. On Windows, installing PyTorch into a `.venv` this deep in a directory tree can exceed the 260-character path limit; enable long paths or create the environment in a shorter directory.

## Running a reconstruction

```powershell
python -m fpm data\usaf\edge_upper_right --method mcfpm
```

or from Python:

```python
from fpm import FPM, output_dir

model = FPM.from_directory("data/usaf/edge_upper_right", method="mcfpm")
result = model.run()
result.save(output_dir("data/usaf/edge_upper_right", "mcfpm"), model=model, summary=True)
```

or step by step in [notebooks/demo_exp.ipynb](notebooks/demo_exp.ipynb), which also shows the raw frames, compares the result with conventional FPM and can run the whole collection.

`save` writes `reconstruction.npz` (object, pupil, LED window centres, recovered array pose, loss traces), `metadata.json` (the files and settings that produced it) and `summary.png`. Results go to `notebooks/output/<dataset>/<method>/`, never under `data/`.

To reconstruct every dataset and compare the results:

```powershell
python tools\run_all.py --methods pie mcfpm      # log and outcomes in notebooks\output\run_all\
python tools\summarize_results.py                # notebooks\output\comparison_pie_mcfpm.png
```


## Repository layout

```text
src/fpm/            Python package: one class, plus stateless mathematics
data/<sample>/<segment>/
                    one directory per acquisition: raw/ + config.yaml + optional config_<method>.yaml
notebooks/          demo_exp.ipynb for any dataset and method; results in notebooks/output/
tools/              run_all.py (every dataset) and summarize_results.py (side-by-side comparison)
tests/              pytest suite, including a line-by-line check against the MATLAB loop
docs/               project page and notes
```

| module | holds |
| --- | --- |
| `model.py` | `FPM`: optics, LEDs, measurements, the ptychographic engine, `pie` and `sa` |
| `solvers.py` | `pcfpm` and `mcfpm`, and the rigid array fit |
| `illumination.py` | LED grid, spiral update order, misaligned wavevectors |
| `optimizers.py` | simulated annealing with the defaults of MATLAB's `simulannealbnd` |
| `simulation.py` | acquisitions of a test object with a known misalignment |
| `config.py`, `data.py` | reading the YAML files and the measurement stacks |
| `fourier.py` | centred FFTs, MATLAB-compatible rounding |
| `visualize.py` | the summary figure |

## Datasets and configuration

`config.yaml` holds facts about the measurement and the instrument: data file, optics, LED array and the position of the patch on the sensor. Solver settings live in `config_<method>.yaml`; every key is optional and defaults to the MATLAB value, listed with an explanation in [data/_template](data/_template). A method file can also override the patch position, which is how `data/root/edge_upper_left` reproduces the 2.28x offset that `FPM_nosa_sa_or_pcsa.m` uses without changing the dataset itself.

Directories are named after where the segment lies on the sensor, as seen in the raw frames (row 0 at the top); the positions were measured by locating each centre frame in the full sensor image.

| dataset | file | segment centre from the sensor centre | used by |
| --- | --- | --- | --- |
| `usaf/center` | `20170604-3.mat` | 0 | paper Fig. 6 |
| `usaf/edge_upper_right` | `20170604-1-400.mat` | 550 px up, 750 px right | paper Fig. 6, `main.m` |
| `root/center` | `20170531-3.mat` (200 x 200 px) | 0 | |
| `root/inner_upper_left` | `20170531-3-300200.mat` | 300 px up, 200 px left | |
| `root/inner_upper_left_far` | `20170531-3-300.mat` | 300 px up, 300 px left | paper Fig. 8 |
| `root/middle_upper_left` | `20170531-3-500400.mat` | 500 px up, 400 px left | |
| `root/edge_upper_left` | `20170531-3-700600.mat` | 700 px up, 600 px left | `FPM_nosa_sa_or_pcsa.m` |
| `root/edge_upper_right` | `20170531-3-550-750.mat` | 550 px up, 750 px right | |

## Tests

```powershell
python -m pytest tests
```

The suite takes about a minute on a CPU. It checks the engine against a NumPy transcription of the MATLAB loop, and every solver on simulated data with a known LED misalignment.

## Reference

```
@Article{Zhou2018OE,
title = {Fast and robust misalignment correction of {F}ourier ptychographic microscopy for full field of view reconstruction},
author = {Ao Zhou and Wei Wang and Ni Chen* and Edmund Y. Lam and Byoungho Lee and Guohai Situ},
journal = {Optics Express},
year = {2018},
volume = {26},
number = {18},
pages = {23661-23674},
month = {Aug.},
doi = {10.1364/oe.26.023661},
}
```

## Contact

Ni Chen (nichen.optics@gmail.com)
