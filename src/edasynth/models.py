"""Trained models, out-of-participant evaluation and metrics (docs/METHOD.md section 6).

All evaluation is leave-one-participant-out (LOSO): the test participant never
contributes to scaling, hyper-parameter choice or fitting. Metrics are computed on
the pooled out-of-fold predictions. See docs/METHOD.md section 6.
"""

from __future__ import annotations

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def makeModel(name: str, cfg: dict, C: float | None = None, seed: int = 0):
    m = cfg["models"]
    if name == "logreg":
        lr = m["logreg"]
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=C if C is not None else lr["CDefault"],
                class_weight=lr["classWeight"],
                max_iter=lr["maxIter"],
                solver=lr.get("solver", "lbfgs"),
            ),
        )
    if name == "rf":
        rf = m["rf"]
        return RandomForestClassifier(
            n_estimators=rf["nEstimators"],
            class_weight=rf["classWeight"],
            random_state=seed,
            n_jobs=rf["nJobs"],
        )
    raise ValueError(name)


def chooseC(X, y, groups, cfg, seed):
    lr = cfg["models"]["logreg"]
    grid = lr["CGrid"]
    if lr["CSelection"] == "fixed" or len(grid) == 1:
        return lr["CDefault"]
    k = lr["innerFolds"]
    ug = np.unique(groups)
    splitter = LeaveOneGroupOut() if k == "loso" else GroupKFold(n_splits=min(int(k), ug.size))
    best, bestC = -np.inf, grid[0]
    for c in grid:
        preds = np.zeros(y.size)
        for tr, te in splitter.split(X, y, groups):
            if y[tr].sum() == 0:
                preds[te] = 0.0
                continue
            mdl = makeModel("logreg", cfg, C=c, seed=seed).fit(X[tr], y[tr])
            preds[te] = mdl.predict_proba(X[te])[:, 1]
        score = average_precision_score(y, preds) if y.sum() else 0.0
        if score > best:
            best, bestC = score, c
    return bestC


def losoPredict(X: np.ndarray, y: np.ndarray, groups: np.ndarray, model: str, cfg: dict, seed: int = 0):
    """Out-of-participant predicted probabilities for every window."""
    preds = np.zeros(y.size)
    chosen = []
    for tr, te in LeaveOneGroupOut().split(X, y, groups):
        if y[tr].sum() == 0:
            preds[te] = 0.0
            continue
        if model == "logreg":
            C = chooseC(X[tr], y[tr], groups[tr], cfg, seed)
            chosen.append(C)
            mdl = makeModel("logreg", cfg, C=C, seed=seed)
        else:
            mdl = makeModel(model, cfg, seed=seed)
        mdl.fit(X[tr], y[tr])
        preds[te] = mdl.predict_proba(X[te])[:, 1]
    return preds, chosen


# ---------------------------------------------------------------- metrics
def auprc(y, s) -> float:
    return float(average_precision_score(y, s)) if y.sum() else float("nan")


def auroc(y, s) -> float:
    return float(roc_auc_score(y, s)) if 0 < y.sum() < y.size else float("nan")


def thresholdAtSpecificity(y, s, spec: float = 0.95) -> float:
    """Score threshold giving the requested specificity on the negatives."""
    neg = np.sort(s[y == 0])
    return float(np.quantile(neg, spec)) if neg.size else float("inf")


def sensitivityAtSpecificity(y, s, spec: float = 0.95) -> float:
    thr = thresholdAtSpecificity(y, s, spec)
    pos = s[y == 1]
    return float((pos > thr).mean()) if pos.size else float("nan")


def bootstrapParticipants(y, sA, sB, groups, nBoot: int, seed: int, metric=auprc):
    """Percentile CI of metric(sB) - metric(sA), resampling participants."""
    rng = np.random.default_rng(seed)
    ug = np.unique(groups)
    idxByG = {g: np.where(groups == g)[0] for g in ug}
    diffs = []
    for iteration in range(nBoot):
        pick = rng.choice(ug, ug.size, replace=True)
        idx = np.concatenate([idxByG[g] for g in pick])
        if y[idx].sum() == 0:
            continue
        diffs.append(metric(y[idx], sB[idx]) - metric(y[idx], sA[idx]))
    d = np.asarray(diffs)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5)), d
