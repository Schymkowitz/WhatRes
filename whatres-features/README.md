# whatres-features

C++ source of `whatres-features`, the command-line program that computes the geometric
features behind WhatRes (apparent resolution in Å and realism 0 to 100, from
coordinates alone). It holds the four commands the WhatRes feature pipeline
calls, and nothing else.

| command | what it computes |
|---|---|
| `pdb_quality_check` | per-residue backbone dihedrals, side-chain and backbone burial, bond lengths and angles, clashes, surface patches and cavities, each scored as a percentile against a reference database |
| `pdb_surface_analysis` | per-chain surface statistics: hydrophobic patches, charge clusters, buried cavities and inter-chain interface voids |
| `PDB_dist` | atom-type-pair distance histograms (60 bins, 0 to 9 Å) and per-residue phi / psi / omega |
| `structure_survey` | builds the binary reference database that `pdb_quality_check` scores against, from a set of experimental structures |

## Build

Needs a C++14 compiler and `make`. The atom parameters come as a precompiled
library, `lib/<platform>/libatomparams.a` (see Citation); `make` stops with a
clear message on a platform that has none.

```bash
make
```

This compiles every `.cpp` in the directory into `./whatres-features` (clang++ on macOS,
g++ on Linux; Linux links `-lpthread`). `make debug` builds `./whatres-features_debug`,
`make install PREFIX=/some/bin` copies the binary, `make clean` removes the
build. On macOS `make ARCH=arm64`, `ARCH=x86_64` or `ARCH=universal` selects
the architecture.

## Run

```bash
whatres-features -do <command> -in <pdb_list.txt> -opt <Options.txt> [-pre <prefix>]
whatres-features -do <command> -help
whatres-features -help
```

- `-in` is a text file with one PDB path per line.
- `-opt` is the options file; `Options_example.txt` lists every option these
  commands read. It starts with `<TITLE>WHATRES_optionfile;` (the older title
  `AGADIR_optionfile` is also accepted). Options left out keep their defaults.
- Output files are written to the current directory, named with `<prefix>`
  (from `-pre` or the `<prefix>` option, default `PSX`).
- Use space-free paths for the list, the PDB files and the working directory:
  the list reader splits on spaces.

`pdb_quality_check` needs a reference database (`<qc_reference_db>`). WhatRes
uses `native49_structure_reference.bin`, part of the WhatRes data deposit. To
build one from scratch:

```bash
whatres-features -do structure_survey -in reference_pdbs.txt -opt Options.txt
```

This writes `<prefix>structure_reference.bin` plus per-residue, cavity, patch
and summary tables.

## Validation

The output of this source tree is byte-identical to the build used for the
WhatRes manuscript: on 12 reference structures, the 608-feature matrix and
every output file of the four commands; and a rebuild of the reference
database from its 2,000 structures.

## Citation

The burial terms use the FoldX occupancy method (`FxAtoms.h`) with the FoldX
atom parameters: per-atom volume and minimal and maximal occupancy. The
parameter values are distributed only in compiled form, under the FoldX
licence (`lib/<platform>/libatomparams.a`, interface in `atomparams.h`).
Please cite FoldX alongside WhatRes:

Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. The FoldX web
server: an online force field. *Nucleic Acids Res.* 2005;33:W382-W388.
doi:10.1093/nar/gki387

## Licence

The source code is released under the BSD 3-Clause licence (`LICENSE`),
copyright VIB. The compiled atom-parameter library is not covered by it; see
`lib/NOTICE.txt`.

## Source

Extracted from the Switch Lab structure toolkit (commit 26f03bb).
