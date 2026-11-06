#!/usr/bin/env python3
"""
build_feature_matrix.py
=======================

Unified per-structure feature extractor for the structure-quality project.

Harvests *everything* whatres-features can currently compute about a
set of PDB structures and merges it into one wide per-structure table
(`feature_matrix.tsv`, one row per structure). This is "Track A" of the
feature-expansion plan: it uses the EXISTING commands, no C++ schema change.

Feature families produced (prefix in column name):
  qc_*      per-residue quality metrics from `pdb_quality_check`
            (36 percentile metrics + 35 surprise scores), aggregated per
            structure as mean / sd / frac_extreme, plus residue-score and
            outlier-rate summaries and secondary-structure composition.
  surf_*    surface patches, charge clusters, cavities, interface voids from
            `pdb_surface_analysis` (aggregated across chains).
  dist_*    pairwise atom-type distance distribution summaries from `PDB_dist`
            (per-structure short/mid/long-range contact fractions per atom-type
            pair, and optional divergence vs a reference histogram).
  bf_*      per-atom B-factor / pLDDT distribution stats parsed straight from
            the PDB coordinates.

Usage
-----
  python3 build_feature_matrix.py \
      --bin   /path/to/whatres-features \
      --list  benchmark/list_all_benchmark.txt \
      --opt   analysis/Options_features_native.txt \
      --root  /path/to/HighResPDB \
      --out   analysis/feature_matrix.tsv \
      --workdir analysis/work

  # fast iteration without the quality check (skips the qc_* family):
  python3 build_feature_matrix.py ... --skip-qc

Notes
-----
* The whatres-features commands take a *list* of PDBs and write combined outputs, so
  `pdb_quality_check` and `pdb_surface_analysis` are each run ONCE on the whole
  list. `PDB_dist` aggregates all inputs
  into a single histogram set, so to get per-structure distance features it is
  run once per structure (cheap).
* Structure id = basename of the PDB path with the extension stripped,
  lowercased. This matches how the commands label rows.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd


# ──────────────────────────────────────────────────────────────────────────────
# helpers
# ──────────────────────────────────────────────────────────────────────────────

def sid(path):
    """Canonical structure id from a path: basename, no extension, lowercase."""
    return os.path.splitext(os.path.basename(path))[0].lower()


def run_cmd(binpath, do, list_file, opt_file, cwd):
    """Run `whatres-features -do <do> -in <list_file> -opt <opt_file>` in cwd."""
    cmd = [binpath, "-do", do, "-in", list_file, "-opt", opt_file]
    res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True)
    # A binary that rejects the options file (missing, or a <TITLE> it does not
    # know) silently runs on its built-in defaults. Stop instead.
    for msg in ("Not a valid optionfile", "Optionfile not found"):
        if msg in res.stdout:
            sys.exit(f"[features] {do}: {binpath} did not accept the options file {opt_file} "
                     f"({msg}); it would have run on default settings. Check the binary and the "
                     f"<TITLE> line of the options file.")
    return res.returncode, res.stdout


def read_list(list_file, root):
    """Read a PDB list file; return list of (id, absolute_path)."""
    out = []
    with open(list_file) as f:
        for line in f:
            p = line.strip()
            if not p or p.startswith("#"):
                continue
            ap = p if os.path.isabs(p) else os.path.join(root, p)
            out.append((sid(p), ap))
    return out


def parse_opt_value(opt_file, tag):
    """Return the value of <tag>value; from an options file, else None."""
    pat = re.compile(r"<" + re.escape(tag) + r">([^;]*);")
    with open(opt_file) as f:
        for ln in f:
            m = pat.search(ln)
            if m:
                return m.group(1).strip()
    return None


def stage_scratch(scratch, structs, opt_file):
    """
    Stage everything under a SPACE-FREE scratch dir via symlinks, because the
    whatres-features PDB-list reader breaks on paths
    containing spaces. Returns (patched_opt, list of
    (id, scratch_pdb_path)).
    """
    if " " in scratch:
        sys.exit(f"--scratch must be a space-free path (got: {scratch!r})")
    pdbs = os.path.join(scratch, "pdbs")
    binl = os.path.join(scratch, "bin")
    os.makedirs(pdbs, exist_ok=True)
    os.makedirs(binl, exist_ok=True)

    staged = []
    seen = set()
    for pid, apath in structs:
        name = pid
        k = name
        n = 1
        while k in seen:                       # de-dup identical ids
            k = f"{name}_{n}"; n += 1
        seen.add(k)
        link = os.path.join(pdbs, k + ".pdb")
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(apath), link)
        staged.append((k, link))

    # symlink the reference db into the space-free location
    patched_lines = []
    ref = parse_opt_value(opt_file, "qc_reference_db")
    if ref and not os.path.isabs(ref):          # relative paths are relative to the options file
        ref = os.path.join(os.path.dirname(os.path.abspath(opt_file)), ref)
    ref_link = None
    if ref and os.path.exists(ref):
        ref_link = os.path.join(binl, "ref.bin")
        if os.path.islink(ref_link) or os.path.exists(ref_link):
            os.remove(ref_link)
        os.symlink(os.path.abspath(ref), ref_link)

    with open(opt_file) as f:
        for ln in f:
            if ref_link and "<qc_reference_db>" in ln:
                ln = f"<qc_reference_db>{ref_link};\n"
            patched_lines.append(ln)
    patched_opt = os.path.join(scratch, "Options.txt")
    with open(patched_opt, "w") as f:
        f.writelines(patched_lines)
    return patched_opt, staged


# ──────────────────────────────────────────────────────────────────────────────
# family: qc_*  (per-residue quality metrics, aggregated)
# ──────────────────────────────────────────────────────────────────────────────

def aggregate_quality_summary(summary_path):
    """Per-structure whole-structure metrics from quality_summary.txt → ws_*
    (incl. ws_clashscore = clashes/1000 atoms, ws_res_score, ws_combined_score)."""
    import io
    if not os.path.exists(summary_path):
        return pd.DataFrame()
    lines = open(summary_path).read().splitlines()
    hdr = next((i for i, l in enumerate(lines) if l.startswith("PDB\t")), None)
    if hdr is None:
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO("\n".join(lines[hdr:])), sep="\t")
    df = df[df["PDB"].notna()].copy()
    df = df[~df["PDB"].astype(str).str.contains(" |---", regex=True)]  # drop summary footer lines
    df["id"] = df["PDB"].astype(str).str.lower().map(lambda s: os.path.splitext(s)[0])
    rename = {"combined_score": "ws_combined_score", "res_score": "ws_res_score",
              "ws_score": "ws_ws_score", "mean_res_score": "ws_mean_res_score",
              "worst5_mean": "ws_worst5_mean", "n_residues": "ws_n_residues",
              "res_outliers": "ws_res_outliers", "clashscore": "ws_clashscore",
              "n_cavities": "ws_n_cavities", "total_cavity_vol": "ws_total_cavity_vol",
              "void_fraction": "ws_void_fraction"}
    have = {k: v for k, v in rename.items() if k in df.columns}
    out = df[["id"] + list(have.keys())].rename(columns=have)
    return out.drop_duplicates("id")


LEGACY_ENERGY_METRICS = {
    "total_energy", "backbone_hbond", "sidechain_hbond", "vdw", "electrostatics", "solv_polar",
    "solv_hydrophobic", "vdw_clashes", "entropy_sc", "entropy_mc", "sloop_entropy", "mloop_entropy",
    "cis_bond", "torsional_clash", "backbone_clash", "helix_dipole", "water_bridge", "disulfide",
    "electrostatic_kon", "partial_covalent", "energy_ionisation", "entropy_complex",
    "hetero_back_hbond", "hetero_side_hbond"}


def aggregate_quality_report(report_tsv):
    """
    Aggregate the per-residue quality_report.tsv to per-structure features.

    For each *_pctile column: mean, sd, frac_extreme (pctile<5 or >95).
    For each *_score  column: mean.
    Plus residue_score mean, n_residues, outlier rate, SS composition.
    """
    if not os.path.exists(report_tsv):
        return pd.DataFrame()
    return aggregate_quality_df(pd.read_csv(report_tsv, sep="\t", low_memory=False))


def aggregate_quality_df(df):
    """aggregate_quality_report on an in-memory per-residue table (e.g. a residue subset)."""
    df = df.copy()
    df["__id"] = df["pdb"].astype(str).str.lower().map(
        lambda s: os.path.splitext(s)[0])

    # pdb_quality_check still writes the 24 legacy energy-metric columns, always at
    # pctile 50 / score 0 (it computes no energies); they carry no information and are
    # never turned into features. Secondary structure is not assigned ('-'), so no
    # SS-composition features either.
    pct_cols = [c for c in df.columns if c.endswith("_pctile")
                and c[:-len("_pctile")] not in LEGACY_ENERGY_METRICS]
    sco_cols = [c for c in df.columns if c.endswith("_score")
                and c not in ("worst_score", "residue_score")
                and c[:-len("_score")] not in LEGACY_ENERGY_METRICS]

    rows = []
    for pid, g in df.groupby("__id"):
        feat = {"id": pid, "qc_n_residues": len(g)}
        # per-residue summary scores
        feat["qc_residue_score_mean"] = g["residue_score"].mean()
        feat["qc_residue_score_sd"] = g["residue_score"].std()
        feat["qc_worst_score_mean"] = g["worst_score"].mean()
        feat["qc_n_outliers_mean"] = g["n_outliers"].mean()
        feat["qc_outlier_rate"] = (g["n_outliers"] > 0).mean()
        # burial distribution
        if "facc" in g:
            feat["qc_facc_mean"] = g["facc"].mean()
            feat["qc_facc_buried_frac"] = (g["facc"] < 0.05).mean()
        # percentile metrics: mean / sd / frac_extreme
        for c in pct_cols:
            v = pd.to_numeric(g[c], errors="coerce")
            base = c[:-len("_pctile")]
            feat[f"qc_{base}_mean_pctile"] = v.mean()
            feat[f"qc_{base}_sd_pctile"] = v.std()
            feat[f"qc_{base}_frac_extreme"] = ((v < 5) | (v > 95)).mean()
        # score metrics: mean
        for c in sco_cols:
            v = pd.to_numeric(g[c], errors="coerce")
            base = c[:-len("_score")]
            feat[f"qc_{base}_mean_score"] = v.mean()
        rows.append(feat)
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# family: surf_*  (surface patches / cavities / interface voids)
# ──────────────────────────────────────────────────────────────────────────────

def aggregate_surface(surf_tsv, cav_tsv, iface_tsv):
    """Aggregate pdb_surface_analysis outputs to per-structure features."""
    feats = defaultdict(dict)

    # per-chain surface analysis -> aggregate across chains
    if os.path.exists(surf_tsv):
        s = pd.read_csv(surf_tsv, sep="\t", low_memory=False)
        s["__id"] = s["pdb_name"].astype(str).str.lower().map(
            lambda x: os.path.splitext(x)[0])
        # sum-like (counts), max-like (largest patch quantities), mean-like
        sum_cols = ["n_res", "n_surf", "n_hphob_patches",
                    "n_pos_clusters", "n_neg_clusters"]
        max_cols = ["hphob1_size", "hphob1_kd_sum", "hphob1_rg",
                    "hphob1_aromatic", "hphob1_entropy",
                    "pos1_size", "pos1_charge", "pos1_rg",
                    "neg1_size", "neg1_charge", "neg1_rg"]
        mean_cols = ["frac_surf", "mean_burial", "surf_net_charge",
                     "frac_surf_hphob", "hphob1_mean_burial",
                     "pos1_mean_burial", "neg1_mean_burial"]
        for pid, g in s.groupby("__id"):
            d = feats[pid]
            d["surf_n_chains"] = len(g)
            for c in sum_cols:
                if c in g:
                    d[f"surf_{c}_sum"] = pd.to_numeric(g[c], errors="coerce").sum()
            for c in max_cols:
                if c in g:
                    d[f"surf_{c}_max"] = pd.to_numeric(g[c], errors="coerce").max()
            for c in mean_cols:
                if c in g:
                    d[f"surf_{c}_mean"] = pd.to_numeric(g[c], errors="coerce").mean()

    # per-pdb cavities
    if os.path.exists(cav_tsv):
        c = pd.read_csv(cav_tsv, sep="\t", low_memory=False)
        c["__id"] = c["pdb_name"].astype(str).str.lower().map(
            lambda x: os.path.splitext(x)[0])
        for pid, g in c.groupby("__id"):
            r = g.iloc[0]
            d = feats[pid]
            for col in ["n_cavities", "total_cavity_vol", "largest_cavity_vol",
                        "protein_volume", "void_fraction"]:
                if col in g:
                    d[f"surf_{col}"] = pd.to_numeric(pd.Series([r[col]]),
                                                     errors="coerce").iloc[0]

    # interface voids (multi-chain only) -> aggregate across chain pairs
    if os.path.exists(iface_tsv):
        try:
            iv = pd.read_csv(iface_tsv, sep="\t", low_memory=False)
        except Exception:
            iv = pd.DataFrame()
        if len(iv) and "pdb_name" in iv:
            iv["__id"] = iv["pdb_name"].astype(str).str.lower().map(
                lambda x: os.path.splitext(x)[0])
            for pid, g in iv.groupby("__id"):
                d = feats[pid]
                d["surf_n_iface_pairs"] = len(g)
                for col in ["iface_void_vol", "packing_density", "void_index"]:
                    if col in g:
                        v = pd.to_numeric(g[col], errors="coerce")
                        d[f"surf_{col}_mean"] = v.mean()

    if not feats:
        return pd.DataFrame()
    out = pd.DataFrame.from_dict(feats, orient="index")
    out.index.name = "id"
    return out.reset_index()


# ──────────────────────────────────────────────────────────────────────────────
# family: dist_*  (atom-type pairwise distance distributions)
# ──────────────────────────────────────────────────────────────────────────────

# distance bin centres in PSX_PDB_dist.txt: 0.075 .. 8.925, width 0.15 (60 bins)
DIST_BIN_W = 0.15
SHORT_HI = 4.0     # contact / short-range
MID_HI = 6.0       # medium-range


def parse_pdb_dist_file(path):
    """Return dict: pair -> (np.array of 60 normalised freqs, total_count)."""
    pairs = {}
    if not os.path.exists(path):
        return pairs
    with open(path) as f:
        header = f.readline()  # 'dist' + bin centres
        centres = np.array([float(x) for x in header.split("\t")[1:] if x.strip()])
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            # last non-empty field is total count; middle 60 are freqs
            vals = [p for p in parts[1:] if p != ""]
            if len(vals) < 2:
                continue
            total = float(vals[-1])
            freqs = np.array([float(v) for v in vals[:-1]])
            pairs[name] = (freqs, total, centres if len(centres) else None)
    return pairs


def dist_features_for_structure(pairs, n_res):
    """
    Summarise one structure's atom-type-pair distance histograms.

    For each pair: total contacts per residue, and fraction of contacts that
    are short (<4A), mid (4-6A), long (6-9A). Column names use the pair label.
    """
    feat = {}
    for name, (freqs, total, centres) in pairs.items():
        if centres is None or freqs.size != centres.size or freqs.sum() <= 0:
            continue
        p = freqs / freqs.sum()
        short = p[centres < SHORT_HI].sum()
        mid = p[(centres >= SHORT_HI) & (centres < MID_HI)].sum()
        long = p[centres >= MID_HI].sum()
        key = name.replace(".", "").replace("_", "-")
        feat[f"dist_{key}_short"] = short
        feat[f"dist_{key}_mid"] = mid
        feat[f"dist_{key}_percontact"] = total / max(n_res, 1)
    return feat


def dist_emd_vs_reference(pairs, ref_pairs):
    """Per-pair 1D earth-mover distance vs a reference histogram + summary."""
    feat = {}
    emds = []
    weights = []
    for name, (freqs, total, centres) in pairs.items():
        if name not in ref_pairs:
            continue
        rfreqs, rtotal, rcent = ref_pairs[name]
        if freqs.sum() <= 0 or rfreqs.sum() <= 0 or freqs.size != rfreqs.size:
            continue
        p = freqs / freqs.sum()
        q = rfreqs / rfreqs.sum()
        emd = np.abs(np.cumsum(p) - np.cumsum(q)).sum() * DIST_BIN_W
        key = name.replace(".", "").replace("_", "-")
        feat[f"dist_{key}_emd"] = emd
        emds.append(emd)
        weights.append(total)
    if emds:
        emds = np.array(emds)
        weights = np.array(weights)
        feat["dist_emd_mean"] = emds.mean()
        feat["dist_emd_max"] = emds.max()
        if weights.sum() > 0:
            feat["dist_emd_wmean"] = float((emds * weights).sum() / weights.sum())
    return feat


# ──────────────────────────────────────────────────────────────────────────────
# family: bf_*  (B-factor / pLDDT from coordinates)
# ──────────────────────────────────────────────────────────────────────────────

def bfactor_features(pdb_path):
    """Per-CA B-factor / pLDDT distribution stats parsed directly from the PDB."""
    cas = []
    allb = []
    try:
        with open(pdb_path, errors="ignore") as f:
            for ln in f:
                if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                    continue
                try:
                    b = float(ln[60:66])
                except ValueError:
                    continue
                allb.append(b)
                if ln[12:16].strip() == "CA":
                    cas.append(b)
    except OSError:
        return {}
    if not allb:
        return {}
    ca = np.array(cas) if cas else np.array(allb)
    ab = np.array(allb)
    feat = {
        "bf_ca_mean": ca.mean(), "bf_ca_sd": ca.std(),
        "bf_ca_min": ca.min(), "bf_ca_max": ca.max(),
        "bf_ca_p10": np.percentile(ca, 10), "bf_ca_p90": np.percentile(ca, 90),
        "bf_all_mean": ab.mean(), "bf_all_sd": ab.std(),
    }
    # AlphaFold/Boltz pLDDT lives in 0-100; flag likely-pLDDT
    feat["bf_looks_plddt"] = int(ca.max() <= 100.0 and ca.min() >= 0.0
                                 and ca.mean() > 30.0 and ca.std() < 30.0)
    feat["bf_ca_frac_lt70"] = (ca < 70).mean()
    return feat


# ──────────────────────────────────────────────────────────────────────────────
# main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, help="path to the whatres-features binary")
    ap.add_argument("--list", required=True, help="PDB list file (one path per line)")
    ap.add_argument("--opt", required=True, help="Options file (reference DB and surface tags)")
    ap.add_argument("--root", default=".", help="root dir for relative PDB paths in --list")
    ap.add_argument("--out", required=True, help="output feature_matrix.tsv")
    ap.add_argument("--workdir", default="analysis/work",
                    help="dir for the output (may contain spaces)")
    ap.add_argument("--scratch", default="/tmp/whatres_featbuild",
                    help="SPACE-FREE scratch dir for staging + command runs")
    ap.add_argument("--skip-qc", action="store_true",
                    help="skip pdb_quality_check (no qc_* features; fast)")
    ap.add_argument("--skip-dist", action="store_true",
                    help="skip PDB_dist per-structure features")
    ap.add_argument("--ref-hist", default=None,
                    help="optional reference PSX_PDB_dist.txt for dist_*_emd features")
    args = ap.parse_args()

    binpath = os.path.abspath(args.bin)
    root = os.path.abspath(args.root)
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)

    raw_structs = read_list(os.path.abspath(args.list), root)
    print(f"[features] {len(raw_structs)} structures from {args.list}", flush=True)

    # Stage into a space-free scratch dir (the whatres-features list reader breaks on
    # spaces).  All commands run under this scratch.
    scratch = os.path.abspath(args.scratch)
    os.makedirs(scratch, exist_ok=True)
    opt_abs, structs = stage_scratch(scratch, raw_structs, os.path.abspath(args.opt))
    print(f"[features] staged {len(structs)} structures + options under {scratch}",
          flush=True)

    # absolute-path list file (now space-free) the commands consume
    abs_list = os.path.join(scratch, "_abs_list.txt")
    with open(abs_list, "w") as f:
        for _, ap_ in structs:
            f.write(ap_ + "\n")

    tables = []

    # ── qc_* ────────────────────────────────────────────────────────────────
    if not args.skip_qc:
        qc_dir = os.path.join(scratch, "qc")
        os.makedirs(qc_dir, exist_ok=True)
        print("[features] running pdb_quality_check…", flush=True)
        rc, log = run_cmd(binpath, "pdb_quality_check", abs_list, opt_abs, qc_dir)
        # locate the report by the configured prefix
        rep = None
        for fn in os.listdir(qc_dir):
            if fn.endswith("quality_report.tsv"):
                rep = os.path.join(qc_dir, fn)
        if rep:
            qc = aggregate_quality_report(rep)
            print(f"[features]   qc_* features for {len(qc)} structures "
                  f"({qc.shape[1]-1} cols)", flush=True)
            tables.append(qc)
            # whole-structure metrics from the summary (ws_clashscore etc.)
            summ = None
            for fn in os.listdir(qc_dir):
                if fn.endswith("quality_summary.txt"):
                    summ = os.path.join(qc_dir, fn)
            if summ:
                ws = aggregate_quality_summary(summ)
                if len(ws):
                    print(f"[features]   ws_* features for {len(ws)} structures "
                          f"({ws.shape[1]-1} cols)", flush=True)
                    tables.append(ws)
        else:
            print(f"[features]   WARNING: no quality_report.tsv produced (rc={rc})")
            print(log[-800:])

    # ── surf_* ───────────────────────────────────────────────────────────────
    surf_dir = os.path.join(scratch, "surf")
    os.makedirs(surf_dir, exist_ok=True)
    print("[features] running pdb_surface_analysis…", flush=True)
    rc, log = run_cmd(binpath, "pdb_surface_analysis", abs_list, opt_abs, surf_dir)
    surf_tsv = cav_tsv = iface_tsv = None
    for fn in os.listdir(surf_dir):
        if fn.endswith("_surface_analysis.txt"):
            surf_tsv = os.path.join(surf_dir, fn)
        elif fn.endswith("_cavities.txt"):
            cav_tsv = os.path.join(surf_dir, fn)
        elif fn.endswith("_interface_voids.txt"):
            iface_tsv = os.path.join(surf_dir, fn)
    surf = aggregate_surface(surf_tsv or "", cav_tsv or "", iface_tsv or "")
    if len(surf):
        print(f"[features]   surf_* features for {len(surf)} structures "
              f"({surf.shape[1]-1} cols)", flush=True)
        tables.append(surf)
    else:
        print(f"[features]   WARNING: no surface output (rc={rc})")

    # ── dist_* and bf_* (per structure) ──────────────────────────────────────
    ref_pairs = None
    if args.ref_hist:
        ref_pairs = parse_pdb_dist_file(os.path.abspath(args.ref_hist))
        print(f"[features] reference histogram: {len(ref_pairs)} atom-pairs")

    dist_rows = []
    bf_rows = []
    dist_dir = os.path.join(scratch, "dist")
    os.makedirs(dist_dir, exist_ok=True)
    for i, (pid, apath) in enumerate(structs):
        # bf_* (cheap, always)
        bf = bfactor_features(apath)
        if bf:
            bf["id"] = pid
            bf_rows.append(bf)
        # dist_*
        if not args.skip_dist:
            one_list = os.path.join(dist_dir, "_one.txt")
            with open(one_list, "w") as f:
                f.write(apath + "\n")
            run_cmd(binpath, "PDB_dist", one_list, opt_abs, dist_dir)
            dfile = os.path.join(dist_dir, "PSX_PDB_dist.txt")
            pairs = parse_pdb_dist_file(dfile)
            # n_res for normalisation
            nres = None
            dih = os.path.join(dist_dir, "PSX_PDB_dihedrals.txt")
            if os.path.exists(dih):
                try:
                    nres = sum(1 for _ in open(dih)) - 1
                except OSError:
                    nres = None
            feat = dist_features_for_structure(pairs, nres or 1)
            if ref_pairs:
                feat.update(dist_emd_vs_reference(pairs, ref_pairs))
            if feat:
                feat["id"] = pid
                dist_rows.append(feat)
        if (i + 1) % 50 == 0:
            print(f"[features]   per-structure dist/bf: {i+1}/{len(structs)}",
                  flush=True)

    if bf_rows:
        tables.append(pd.DataFrame(bf_rows))
    if dist_rows:
        d = pd.DataFrame(dist_rows)
        print(f"[features]   dist_* features for {len(d)} structures "
              f"({d.shape[1]-1} cols)", flush=True)
        tables.append(d)

    # ── merge all families on id ─────────────────────────────────────────────
    if not tables:
        sys.exit("[features] no feature tables produced, aborting")
    merged = tables[0]
    for t in tables[1:]:
        merged = merged.merge(t, on="id", how="outer")
    merged = merged.sort_values("id").reset_index(drop=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"[features] wrote {args.out}: {merged.shape[0]} structures × "
          f"{merged.shape[1]-1} features", flush=True)


if __name__ == "__main__":
    main()
