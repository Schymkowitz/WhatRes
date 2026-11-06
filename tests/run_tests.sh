#!/usr/bin/env bash
# Build whatres-features, score the 12 test structures, compare features and scores with the expected values.
set -euo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"; PY="${PYTHON:-python3}"
make -C "$R/whatres-features" >/dev/null 2>&1 || { make -C "$R/whatres-features"; exit 1; }
W="$(mktemp -d /tmp/whatres_test.XXXXXX)"
cd "$R" && "$PY" -m whatres.cli score tests/structures/*.pdb -o "$W/scores.tsv" --features-out "$W/features.tsv" >/dev/null
"$PY" - "$R/tests" "$W" <<'PY'
import sys, numpy as np, pandas as pd
T, W = sys.argv[1:]
def cmp(a, b, key):
    a, b = pd.read_csv(a, sep="\t").set_index(key).sort_index(), pd.read_csv(b, sep="\t").set_index(key).sort_index()
    assert list(a.index) == list(b.index), "different structures"
    num = [c for c in a.columns if c in b.columns and pd.api.types.is_numeric_dtype(a[c])]
    d = max((np.abs(a[c] - b[c]) / np.maximum(np.abs(a[c]), 1e-12)).max() for c in num)
    return len(a), len(num), set(a.columns) ^ set(b.columns), d
n, m, extra, d = cmp(f"{T}/expected_features.tsv", f"{W}/features.tsv", "id")
print(f"features: {n} structures x {m} columns, max relative difference {d:.1e}, column differences {sorted(extra)}")
ok = d < 1e-6 and not extra
n, m, extra, d2 = cmp(f"{T}/expected_scores.tsv", f"{W}/scores.tsv", "id")
print(f"scores:   {n} structures, max relative difference {d2:.1e}")
ok = ok and d2 < 1e-6
print("PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
PY
