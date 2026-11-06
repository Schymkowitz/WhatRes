#!/usr/bin/env python3
"""WhatRes command line.

  python -m whatres.cli score <structure.pdb ...> [-o scores.tsv] [--features-out features.tsv]
  python -m whatres.cli fetch-models

score: computes the features with whatres-features ($WHATRES_FEATURES, else the binary built in
whatres-features/, else on PATH), then apparent resolution and realism. PDB-format input.

Training structures: the apparent-resolution model reproduces the deposited resolution of the
structures it was trained on. For an input whose file name is a training-structure PDB id (for
example 1abc.pdb), the out-of-fold value from the published cross-validation is reported instead,
and the `note` column says so. Renamed copies of training structures are not recognised."""
import argparse, hashlib, os, shutil, subprocess, sys, tempfile, urllib.request
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, MODELS = os.path.join(ROOT, "data"), os.environ.get("WHATRES_MODELS_DIR", os.path.join(ROOT, "models"))
MODEL_URL = "https://github.com/Schymkowitz/WhatRes/releases/download/v1.0.0/apparent_resolution_model.joblib"
MODEL_SHA256 = "a1e5d0ba4433b3202e93486c9519465855ca7a0863efc226854cef27777a1473"

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def fetch_models():
    dst = os.path.join(MODELS, "apparent_resolution_model.joblib")
    if os.path.exists(dst) and sha256(dst) == MODEL_SHA256:
        print(f"model present and verified: {dst}"); return
    print(f"downloading {MODEL_URL}"); tmp = dst + ".part"
    urllib.request.urlretrieve(MODEL_URL, tmp)
    if sha256(tmp) != MODEL_SHA256: os.remove(tmp); sys.exit("checksum mismatch: download rejected")
    os.replace(tmp, dst); print(f"model saved and verified: {dst}")

def binary():
    for b in (os.environ.get("WHATRES_FEATURES"), os.path.join(ROOT, "whatres-features", "whatres-features"), shutil.which("whatres-features")):
        if b and os.path.exists(b): return b
    sys.exit("whatres-features not found: build it (make -C whatres-features) or set $WHATRES_FEATURES")

def score(files, out, features_out):
    if not os.path.exists(os.path.join(MODELS, "apparent_resolution_model.joblib")):
        sys.exit("apparent-resolution model missing: run `python -m whatres.cli fetch-models`")
    from . import features, scorer
    with tempfile.TemporaryDirectory(prefix="whatres_") as td:          # space-free scratch
        opt = os.path.join(td, "Options.txt"); ref = os.path.join(DATA, "native49_structure_reference.bin")
        s = open(os.path.join(DATA, "Options_features.txt")).read()
        s = "\n".join(f"<qc_reference_db>{ref};" if l.startswith("<qc_reference_db>") else l for l in s.split("\n"))
        open(opt, "w").write(s)
        lst = os.path.join(td, "list.txt"); open(lst, "w").write("\n".join(os.path.abspath(f) for f in files) + "\n")
        fm = features_out or os.path.join(td, "features.tsv")
        sys.argv = ["features", "--bin", binary(), "--list", lst, "--opt", opt, "--root", "/", "--out", fm,
                    "--workdir", os.path.join(td, "work"), "--scratch", os.path.join(td, "scratch")]
        features.main()
        df = pd.read_csv(fm, sep="\t")
    ids = df["id"].astype(str).str.lower().str.replace(r"\.(pdb|ent)$", "", regex=True)
    ap = scorer.apparent_resolution(df); rl = scorer.realism_score(scorer.realism_load(df))
    oa = pd.read_csv(os.path.join(MODELS, "training_oof_apparent_resolution.tsv"), sep="\t").set_index("id")["oof"]
    orl = pd.read_csv(os.path.join(MODELS, "training_oof_realism.tsv"), sep="\t").set_index("id")["realism_oof"]
    res = pd.DataFrame({"id": ids, "apparent_resolution": ap, "realism": rl, "note": ""})
    t = res.id.isin(oa.index) | res.id.isin(orl.index)
    res.loc[res.id.isin(oa.index), "apparent_resolution"] = res.id.map(oa)
    res.loc[res.id.isin(orl.index), "realism"] = res.id.map(orl)
    res.loc[t, "note"] = "training structure: out-of-fold value"
    res.to_csv(out or sys.stdout, sep="\t", index=False, float_format="%.4f")

def main():
    ap = argparse.ArgumentParser(prog="whatres", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("score"); s.add_argument("structures", nargs="+"); s.add_argument("-o", "--out"); s.add_argument("--features-out")
    sp.add_parser("fetch-models")
    a = ap.parse_args()
    fetch_models() if a.cmd == "fetch-models" else score(a.structures, a.out, a.features_out)

if __name__ == "__main__":
    main()
