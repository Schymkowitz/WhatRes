# WhatRes

WhatRes scores a protein structure from its coordinates alone, on two axes:

- **apparent resolution** (Å): the deposited resolution at which an experimental structure would
  typically show the same local geometry;
- **realism** (0 to 100): how experimental-like the geometry is, independent of resolution
  (100 = experimental-like, low = over-idealised, as in many predicted models).

It works on experimental structures (X-ray, cryo-EM, NMR) and predicted models alike, and needs no
experimental data and no energy function.

## Install

Requirements: a C++14 compiler and `make`; Python 3.11 or later.

```bash
git clone https://github.com/Schymkowitz/WhatRes.git && cd WhatRes
make -C whatres-features                 # builds whatres-features/whatres-features
python3 -m pip install -r requirements.txt
python3 -m whatres.cli fetch-models      # apparent-resolution model (77 MB, checksum-verified)
tests/run_tests.sh                       # 12 test structures: expect PASS
```

`make` links the precompiled atom-parameter library for your platform from `whatres-features/lib/`;
the platforms provided are the subfolders there.

## Use

```bash
python3 -m whatres.cli score my_structure.pdb other.pdb -o scores.tsv
```

Input: PDB-format files (one model is read). Output, one row per structure:

| column | meaning |
|---|---|
| `id` | file name without extension |
| `apparent_resolution` | Å |
| `realism` | 0 to 100 |
| `note` | `training structure: out-of-fold value` when applicable (below) |

`--features-out features.tsv` also writes the per-structure feature table.

From Python, with a feature table:

```python
import pandas as pd, whatres
df = pd.read_csv("features.tsv", sep="\t")
ap = whatres.apparent_resolution(df)
realism = whatres.realism_score(whatres.realism_load(df))
```

## Training structures

The apparent-resolution model reproduces the deposited resolution of the structures it was trained
on, so a score for one of them would not be a prediction. When an input file is named by the PDB id
of a training structure (for example `1abc.pdb`), WhatRes reports the out-of-fold value from the
published cross-validation instead and says so in the `note` column. Copies of training structures
under other names are not recognised.

## Contents

| path | |
|---|---|
| `whatres-features/` | C++ source of the feature program (see its README) |
| `whatres/` | Python: feature assembly (`features.py`), scoring (`scorer.py`), command line (`cli.py`) |
| `data/` | options file and the high-resolution reference distributions |
| `models/` | realism model and the out-of-fold values of the training structures |
| `tests/` | test structures with their expected features and scores |

The code that trains the models and reproduces the analyses of the paper is in
[WhatRes-paper](https://github.com/Schymkowitz/WhatRes-paper).

## Citation

If you use WhatRes, please cite the WhatRes paper (see `CITATION.cff`). The burial terms use the
per-atom volumes and occupancy limits of the FoldX 2014 parameter set, distributed here in compiled
form (`whatres-features/lib/NOTICE.txt`); please also cite:
Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web server: an online
force field. *Nucleic Acids Res.* 2005;33:W382-W388. doi:10.1093/nar/gki387

The all-atom clashscore uses the van der Waals radii of MolProbity's probe (`eRad`, atomprops.h;
https://github.com/rlabduke/probe, Apache 2.0).

## Licence

BSD 3-Clause (`LICENSE`), except the compiled atom-parameter library (`whatres-features/lib/NOTICE.txt`).
