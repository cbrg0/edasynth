"""Validation gates G2-G7 (G1 is the pytest suite).

Thresholds come from ``configs/*.yaml`` (section ``validation``) and are fixed
before the first run. Each gate returns a dict with ``passed`` and the numbers
behind the decision; the report is written by ``cli validate``.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .data import Recording
from .estimation import decomposeParticipant, responseRates, scrEvents
from .experiment import evaluate, makeDataset
from .features import f0
from .generator import Condition, Generator, assignDonors


# ------------------------------------------------------------------ G2
def g2ParameterRecovery(bank: dict, cfg: dict, n: int, seed: int) -> dict:
    """Generate artifact-free participants with known physiology; re-estimate it."""
    v = cfg["validation"]
    rng = np.random.default_rng(seed)
    gen = Generator(bank["library"], cfg)
    cond = Condition(prevalence=0.0, piR=0.0)
    mAll, pAll = bank["motion"], bank["phys"]
    perm = assignDonors(mAll, pAll, rng)
    tauErr, rateRel, trueA, estA = [], [], [], []
    for k in range(min(n, len(mAll))):
        m, p = mAll[k], pAll[int(perm[k])]
        rec = gen.generate(f"g2-{k}", m, p, cond, rng)
        nw = rec.labels.size
        r = Recording(pid=rec.pid, fs=rec.fs, eda=rec.y, acc=rec.acc, stage=rec.stage, labels3=np.zeros((nw, 3), int))
        clean = np.ones(nw, bool)
        dec = decomposeParticipant(r, clean, cfg)
        ev = scrEvents(
            dec["r"], r.fs, dec["tau0"], cfg["estimation"]["cvxeda"]["tau1"], cfg["estimation"]["minResponseAmp"]
        )
        rates = responseRates(r, ev, clean)
        tauErr.append(abs(dec["tau0"] - p["tau0"]))
        for g in ("rest", "load"):
            t = sum(
                1
                for x in rec.responses
                if x["origin"] == "spontaneous" and x["amp"] >= 0.05 and stageTypeAt(rec, x["t"]) == g
            )
            T = rates[g]["TS"]
            if T > 0 and t > 0:
                rateRel.append(rates[g]["ratePerS"] / (t / T) - 1)
        trueA += [x["area"] for x in rec.responses if x["amp"] >= 0.05]
        estA += [e["area"] for e in ev]
    ks = stats.ks_2samp(trueA, estA) if trueA and estA else None
    medTau = float(np.median(tauErr))
    medRate = float(np.median(np.abs(rateRel))) if rateRel else float("nan")
    ok = medTau <= v["g2Tau0MaxErrS"] and medRate <= v["g2RateMaxRel"] and ks is not None and ks.pvalue > v["g2KsP"]
    return {
        "passed": bool(ok),
        "medianTau0ErrorS": medTau,
        "medianAbsRateRelError": medRate,
        "ksAreaP": float(ks.pvalue) if ks else float("nan"),
        "ksAreaD": float(ks.statistic) if ks else float("nan"),
        "nParticipants": len(tauErr),
    }


def stageTypeAt(rec, t):
    from .data import stageType

    i = min(int(t * rec.fs), rec.stage.size - 1)
    return stageType(str(rec.stage[i]))


# ------------------------------------------------------------------ G3, G4
def classifierTwoSample(Xr, gr, Xs, gs, folds: int = 5) -> dict:
    """AUC of a classifier separating real from synthetic windows (grouped CV)."""
    X = np.vstack([Xr, Xs])
    y = np.r_[np.zeros(len(Xr)), np.ones(len(Xs))]
    g = np.r_[gr, gs + 10000]
    pred = np.zeros(y.size)
    for tr, te in GroupKFold(n_splits=folds).split(X, y, g):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    D = [float(stats.ks_2samp(Xr[:, k], Xs[:, k]).statistic) for k in range(X.shape[1])]
    return {"auc": float(roc_auc_score(y, pred)), "ksD": D}


def g3G4Realism(real: list[Recording], bank: dict, cfg: dict, seed: int) -> dict:
    v = cfg["validation"]
    deriv = cfg["features"]["derivative"]
    mix = tuple(cfg["grid"]["compositions"][-1])
    ds = makeDataset(bank, Condition(c=0.5, classes=mix, prevalence=cfg["grid"]["prevalence"]), cfg, seed)
    Xs, ys, gs = ds["F0"], ds["y"], ds["groups"]
    Xr, yr, gr, unan = [], [], [], []
    for k, rec in enumerate(real):
        F = f0(rec.eda, deriv)
        nw = min(F.shape[0], rec.labels3.shape[0])
        Xr.append(F[:nw])
        yr.append(rec.majority()[:nw])
        unan.append(rec.labels3[:nw].sum(1) == 0)
        gr.append(np.full(nw, k))
    Xr, yr, gr, unan = np.vstack(Xr), np.concatenate(yr), np.concatenate(gr), np.concatenate(unan)
    g3 = classifierTwoSample(Xr[unan], gr[unan], Xs[ys == 0], gs[ys == 0])
    g4 = classifierTwoSample(Xr[yr == 1], gr[yr == 1], Xs[ys == 1], gs[ys == 1])
    g3["passed"] = g3["auc"] <= v["c2stAucMax"]
    g4["passed"] = g4["auc"] <= v["c2stAucMax"]
    g4["note"] = "UTD artifact labels partly depend on movement; interpret with that caveat"
    return {"G3": g3, "G4": g4}


# ------------------------------------------------------------------ G5, G6
def runControlReps(bank, cond, cfg, reps, seed0, accShift=False):
    rows = []
    for r in range(reps):
        ds = makeDataset(bank, cond, cfg, seed0 + r, accShift=accShift)
        row = evaluate(ds, cfg, ("logreg",), seed0 + r)[0]
        row.pop("predictions")
        rows.append(row)
    return rows


def tCi(x):
    x = np.asarray(x, float)
    m = x.mean()
    h = stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size)
    return float(m), float(m - h), float(m + h)


def g5NegativeControls(bank, cfg, seed) -> dict:
    v = cfg["validation"]
    mix = tuple(cfg["grid"]["compositions"][-1])
    out = {}
    ok = True
    for name, cond, shift in (
        ("c0", Condition(c=0.0, classes=mix), False),
        ("accShift", Condition(c=1.0, kappa=1.0, deltaMax=1.0, classes=mix), True),
    ):
        rows = runControlReps(bank, cond, cfg, v["controlReps"], seed, accShift=shift)
        res = {}
        for d in ("dA", "dA4", "dB"):
            m, lo, hi = tCi([r[d] for r in rows])
            good = lo <= 0 <= hi and abs(m) < v["g5MaxAbsDelta"]
            res[d] = {"mean": m, "ci": [lo, hi], "passed": bool(good)}
            ok &= good
        out[name] = res
    out["passed"] = bool(ok)
    return out


def g6PositiveControls(bank, cfg, seed) -> dict:
    v = cfg["validation"]
    out = {}
    ok = True
    for name, cond in (
        ("c4Only", Condition(c=1.0, classes=("C4",))),
        (
            "c1Kappa1Dmax1Mixture",
            Condition(c=1.0, kappa=1.0, deltaMax=1.0, classes=tuple(cfg["grid"]["compositions"][-1])),
        ),
    ):
        rows = runControlReps(bank, cond, cfg, v["controlReps"], seed)
        res = {"nPositiveWindows": [r["nPositive"] for r in rows]}
        anyOk = False
        if min(res["nPositiveWindows"]) == 0:
            res["note"] = "no positive windows were generated; the control cannot be evaluated"
        for d in ("dA", "dB"):
            m, lo, hi = tCi([r[d] for r in rows])
            res[d] = {"mean": m, "ci": [lo, hi]}
            anyOk |= m > 0 and lo > 0
        res["passed"] = bool(anyOk)
        ok &= anyOk
        out[name] = res
    out["passed"] = bool(ok)
    return out


def g7RealConsistency(rawResult: dict, cfg: dict) -> dict:
    v = cfg["validation"]
    ref = v["g7Reference"]
    sameSign = np.sign(rawResult["delta"]) == np.sign(ref["delta"])
    close = (
        abs(rawResult["auprc0"] - ref["auprc0"]) <= v["g7MaxAbsDiff"]
        and abs(rawResult["auprc1"] - ref["auprc1"]) <= v["g7MaxAbsDiff"]
    )
    return {
        "passed": bool(sameSign and close),
        "ours": {k: rawResult[k] for k in ("auprc0", "auprc1", "delta")},
        "reference": ref,
    }
