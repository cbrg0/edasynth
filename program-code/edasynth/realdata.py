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

from .models import auprc, auroc, bootstrapParticipants, losoPredict


def compare(X0, X1, y, g, cfg, model="logreg", seed=0) -> dict:
    p0, unused = losoPredict(X0, y, g, model, cfg, seed)
    p1, unused = losoPredict(X1, y, g, model, cfg, seed)
    lo, hi, boot = bootstrapParticipants(y, p0, p1, g, cfg["analysis"]["bootstrapN"], seed)
    per = []
    for k in np.unique(g):
        m = g == k
        if y[m].sum() > 0:
            per.append(auprc(y[m], p1[m]) - auprc(y[m], p0[m]))
    per = np.asarray(per)
    w = stats.wilcoxon(per) if per.size > 5 else None
    return {
        "model": model,
        "auprc0": auprc(y, p0),
        "auprc1": auprc(y, p1),
        "delta": auprc(y, p1) - auprc(y, p0),
        "ci95": [lo, hi],
        "auroc0": auroc(y, p0),
        "auroc1": auroc(y, p1),
        "prevalence": float(y.mean()),
        "nWindows": int(y.size),
        "nParticipants": int(np.unique(g).size),
        "perParticipantDeltaMedian": float(np.median(per)) if per.size else float("nan"),
        "perParticipantN": int(per.size),
        "wilcoxonP": float(w.pvalue) if w is not None else float("nan"),
        "scores0": p0,
        "scores1": p1,
    }


def shuffledControl(xEda, xAcc, y, g, cfg, nShuffles: int, seed: int) -> dict:
    """Permute the acceleration rows within each participant: alignment is destroyed."""
    rng = np.random.default_rng(seed)
    p0, unused = losoPredict(xEda, y, g, "logreg", cfg, seed)
    base = auprc(y, p0)
    deltas = []
    for iteration in range(nShuffles):
        Xs = xAcc.copy()
        for k in np.unique(g):
            idx = np.where(g == k)[0]
            Xs[idx] = xAcc[rng.permutation(idx)]
        p1, unused = losoPredict(np.hstack([xEda, Xs]), y, g, "logreg", cfg, seed)
        deltas.append(auprc(y, p1) - base)
    return {"meanDelta": float(np.mean(deltas)), "sdDelta": float(np.std(deltas)), "n": nShuffles}


def singleFeatureAuprc(X, y) -> list[float]:
    """AUPRC of each feature used alone as a score (larger value = more artifact-like)."""
    return [auprc(y, X[:, k]) for k in range(X.shape[1])]
