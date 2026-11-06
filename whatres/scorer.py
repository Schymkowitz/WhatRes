"""whatres.py: shared WhatRes scorers. Import and call; the saved models load once.

  apparent_resolution(df) : Angstrom; gradient boosting on the resolution-scaling descriptors
                            (no penalty, no defect exclusion).
  realism_load(df)        : over-regularisation probability P(model-like) in [0,1] from the
                            deployed over_regularisation_model (0 = experimental-like; 1 = over-idealised).
  realism_grade(D)        : 1..6 grade (1 = best, experimental-like).
  realism_score(D)        : 0..100 realism (100 = best, experimental-like) = 100*(1-P).

df is a per-structure feature table (columns = descriptor names). Missing columns are treated
as absent (filled with the training default), so any structure set can be scored uniformly.
"""
import os, joblib, numpy as np, pandas as pd

_H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Model directory is configuration: WHATRES_MODELS_DIR overrides the default analysis/models.
MODELS_DIR = os.environ.get("WHATRES_MODELS_DIR", os.path.join(_H, "models"))
_AP = joblib.load(os.path.join(MODELS_DIR, "apparent_resolution_model.joblib"))
_AP["regressor"].n_jobs = 1   # deterministic: parallel summation of the trees varies at 1e-15 and breaks exact ties
# Deployed realism = the over-regularisation model used for the WhatRes manuscript figures:
# realism = 100 * (1 - P(model-like)). This REPLACES the legacy realism_model.joblib
# (100*exp(-violation_load)), which gave a different, softer scale and is no longer used.
_OR = joblib.load(os.path.join(MODELS_DIR, "over_regularisation_model.joblib"))

SCALE_FEATURES    = _AP["scale_features"]                    # resolution-scaling class (apparent resolution)
GEOMETRY_FEATURES = _AP.get("geometry_features", SCALE_FEATURES)  # size-invariant geometry (UMAP / manifold basis)
REALISM_FEATURES  = _OR["panel"]   # residualised geometry panel feeding the realism classifier

def apparent_resolution(df):
    X = df.reindex(columns=SCALE_FEATURES, fill_value=0.0).replace([np.inf, -np.inf], np.nan).fillna(0.0).values
    return _AP["regressor"].predict(X)

def realism_load(df):
    """Over-regularisation probability P(model-like) in [0,1] from the deployed realism model.
    Each geometry descriptor is residualised on apparent resolution, then an L2 logistic
    regression separates experimental structures (P->0) from predicted models (P->1).
    (0 = experimental-like geometry; 1 = strongly over-idealised.)"""
    ap = np.asarray(apparent_resolution(df))
    cols, CO = _OR["panel"], _OR["appres_resid_coef"]
    sub = df.reindex(columns=cols)
    R = np.column_stack([pd.to_numeric(sub[c], errors="coerce").values - resid_expected(CO[c], ap) for c in cols])
    R = np.nan_to_num(R, nan=0.0)
    return _OR["clf"].predict_proba(R)[:, 1]

def realism_grade(D):
    """1..6 grade (1 = best / experimental-like) from the over-regularisation probability D."""
    return 1 + 5 * np.clip(np.asarray(D, dtype=float), 0, 1)

def realism_score(D):
    """0..100 realism (100 = best / experimental-like) from the over-regularisation probability D."""
    return 100 * (1 - np.clip(np.asarray(D, dtype=float), 0, 1))

def resid_expected(spec, x):
    """Expected (resolution-driven) value of a residualised descriptor at apparent resolution x.
    spec is the entry stored in the realism model's appres_resid_coef:
      ("expsat", [a, b, k])  monotone exponential-saturating fit  a + b*exp(-k*x)  (current)
      ("poly",   coef)       polynomial coefficients for np.polyval (legacy)
      bare ndarray/list      legacy polynomial coefficients (back-compat)
    """
    x = np.asarray(x, dtype=float)
    if isinstance(spec, (tuple, list)) and len(spec) == 2 and isinstance(spec[0], str):
        kind, p = spec
        if kind == "expsat":
            a, b, k = p
            return a + b * np.exp(-k * x)
        return np.polyval(p, x)
    return np.polyval(spec, x)
