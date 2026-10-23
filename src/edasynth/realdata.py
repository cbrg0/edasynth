"""Analyses of the real UTD data (no synthetic data involved).

``stored``  uses the feature matrices published by Zhang et al. (z-scored within
            participant, first derivative by forward difference) and reproduces the
            baseline comparison EDA-only (24 features) vs EDA + 4 acceleration-
            magnitude features (28), logistic regression, leave-one-participant-out.
``raw``     recomputes the features from the PhysioNet signals with this
            package (validation check G7).
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .models import auprc, auroc, bootstrap_participants, loso_predict


def compare(X0, X1, y, g, cfg, model="logreg", seed=0) -> dict:
    p0, _ = loso_predict(X0, y, g, model, cfg, seed)
    p1, _ = loso_predict(X1, y, g, model, cfg, seed)
    lo, hi, boot = bootstrap_participants(y, p0, p1, g, cfg["analysis"]["bootstrap_n"], seed)
    per = []
    for k in np.unique(g):
        m = g == k
        if y[m].sum() > 0:
            per.append(auprc(y[m], p1[m]) - auprc(y[m], p0[m]))
    per = np.asarray(per)
    w = stats.wilcoxon(per) if per.size > 5 else None
    return {
        "model": model,
        "auprc_0": auprc(y, p0),
        "auprc_1": auprc(y, p1),
        "delta": auprc(y, p1) - auprc(y, p0),
        "ci95": [lo, hi],
        "auroc_0": auroc(y, p0),
        "auroc_1": auroc(y, p1),
        "prevalence": float(y.mean()),
        "n_windows": int(y.size),
        "n_participants": int(np.unique(g).size),
        "per_participant_delta_median": float(np.median(per)) if per.size else float("nan"),
        "per_participant_n": int(per.size),
        "wilcoxon_p": float(w.pvalue) if w is not None else float("nan"),
        "_p0": p0,
        "_p1": p1,
    }


def shuffled_control(X_eda, X_acc, y, g, cfg, n_shuffles: int, seed: int) -> dict:
    """Permute the acceleration rows within each participant: alignment is destroyed."""
    rng = np.random.default_rng(seed)
    p0, _ = loso_predict(X_eda, y, g, "logreg", cfg, seed)
    base = auprc(y, p0)
    deltas = []
    for _ in range(n_shuffles):
        Xs = X_acc.copy()
        for k in np.unique(g):
            idx = np.where(g == k)[0]
            Xs[idx] = X_acc[rng.permutation(idx)]
        p1, _ = loso_predict(np.hstack([X_eda, Xs]), y, g, "logreg", cfg, seed)
        deltas.append(auprc(y, p1) - base)
    return {"mean_delta": float(np.mean(deltas)), "sd_delta": float(np.std(deltas)), "n": n_shuffles}


def single_feature_auprc(X, y) -> list[float]:
    """AUPRC of each feature used alone as a score (larger value = more artifact-like)."""
    return [auprc(y, X[:, k]) for k in range(X.shape[1])]
