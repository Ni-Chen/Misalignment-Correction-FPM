# Dataset layout

Each acquisition has one directory, grouped by sample and named after where
the segment lies on the sensor (`center`, `edge_upper_left`, ...):

```text
data/<sample>/<segment>/
├── raw/<file>.mat        # the measurement exactly as recorded (imSeqLowRes, H x W x 289)
├── config.yaml           # data layout, optics, LED array, patch position
└── config_<method>.yaml  # optional: pie, sa, pcfpm, mcfpm
```

To add one, copy `_template/config.yaml`, put the measurement under `raw/`
and fill in `data.file` and the patch position. A method file is only needed
to change a solver setting; every key it leaves out takes the default listed
in `_template/config_<method>.yaml`.

Files are never modified, normalised or cropped by the loader. Results are
written to `notebooks/output/`, not here.

`_assets/` holds the two USC-SIPI images (4.2.03 mandrill, 5.2.09 aerial)
from which `fpm.simulation.load_object` builds the simulated test object.
