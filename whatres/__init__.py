"""WhatRes: apparent resolution (Å) and realism (0 to 100) of a protein structure from its coordinates.

  whatres.apparent_resolution(df)   apparent resolution from a per-structure feature table
  whatres.realism_load(df)          P(model-like) in [0, 1]
  whatres.realism_score(P)          realism 0..100 (100 = experimental-like)
  whatres.features                  feature table from structures (runs whatres-features)
  python -m whatres.cli score ...   the command-line tool

The scoring models load on first use, so `whatres.cli fetch-models` and `whatres.features` work
before the apparent-resolution model has been downloaded."""
__version__ = "1.0.0"

def __getattr__(name):            # lazy: load the models only when a scorer is first used
    import importlib
    if name in ("scorer", "features", "cli"): return importlib.import_module(f".{name}", __name__)
    if name.startswith("__"): raise AttributeError(name)
    sc = importlib.import_module(".scorer", __name__)
    if hasattr(sc, name): return getattr(sc, name)
    raise AttributeError(f"module 'whatres' has no attribute {name!r}")
