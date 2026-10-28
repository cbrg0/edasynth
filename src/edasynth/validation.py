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
    """G2: recovery with cvxEDA's response shape and with varied response shapes; both must pass."""
    plain = g2Recovery(bank, cfg, n, seed, None)
    jitter = g2Recovery(bank, cfg, n, seed + 1, cfg["validation"]["g2Jitter"])
    return {**plain, "variedShapes": jitter, "passed": bool(plain["passed"] and jitter["passed"])}


def g2Recovery(bank: dict, cfg: dict, n: int, seed: int, jitter) -> dict:
    """Generate artifact-free participants with known physiology; re-estimate it.

    With ``jitter`` every synthetic response gets its own time constants, so the data
    no longer have the response shape that the decomposition assumes (docs/METHOD.md
    section 7).
    """
    v = cfg["validation"]
    rng = np.random.default_rng(seed)
    gen = Generator(bank["library"], cfg)
    gen.responseJitter = jitter
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
    """G3/G4 realism without letting the generator see the real windows it is judged on.

    ``validation.realism: loso`` (decision, 1.4.0): for every real participant, one
    synthetic participant is built from the other 19 only (donors and artifact
    episodes); real and synthetic windows of all participants are then pooled in one
    classifier two-sample test. The split-half comparison (generator from one half,
    real windows of the other, both ways) is reported as well. ``realism: pooled``
    is the comparison of version 1.3.0, in which the synthetic artifacts came from
    the same episodes they were compared with.
    """
    scheme = cfg["validation"].get("realism", "pooled")
    if scheme == "pooled":
        return realismPair(real, bank, cfg, seed)
    if scheme != "loso":
        raise ValueError(f"unknown realism scheme {scheme}")
    v = cfg["validation"]
    reps = [realismLoso(real, bank, cfg, seed + 1000 * r) for r in range(v.get("realismReps", 1))]
    out = {}
    for g in ("G3", "G4"):
        aucs = [r[g]["auc"] for r in reps]
        out[g] = {
            "auc": float(np.median(aucs)),
            "aucReplicates": aucs,
            "ksD": list(np.median([r[g]["ksD"] for r in reps], axis=0)),
            "scheme": "loso",
        }
        out[g]["passed"] = bool(out[g]["auc"] <= v["c2stAucMax"])
    split = realismSplitHalf(real, bank, cfg, seed)
    for g in ("G3", "G4"):
        out[g]["splitHalf"] = split[g]
    out["G4"]["note"] = "UTD artifact labels partly depend on movement; interpret with that caveat"
    return out


def realismLoso(real: list[Recording], bank: dict, cfg: dict, seed: int) -> dict:
    """One leave-one-participant-out realism comparison (one synthetic participant per real one)."""
    deriv = cfg["features"]["derivative"]
    mix = tuple(cfg["grid"]["compositions"][-1])
    parts = {k: [] for k in ("Xr3", "gr3", "Xr4", "gr4", "Xs3", "gs3", "Xs4", "gs4")}
    for k, r in enumerate(sorted(real, key=lambda x: x.pid)):
        sub = withoutParticipants(bank, {r.pid}, cfg)
        ds = makeDataset(sub, Condition(c=0.5, classes=mix, prevalence=cfg["grid"]["prevalence"]), cfg, seed + k, 1)
        F = f0(r.eda, deriv)
        nw = min(F.shape[0], r.labels3.shape[0])
        art = r.majority()[:nw] == 1
        clean = r.labels3[:nw].sum(1) == 0
        for key, X, m in (("3", F[:nw], clean), ("4", F[:nw], art)):
            parts["Xr" + key].append(X[m])
            parts["gr" + key].append(np.full(int(m.sum()), k))
        for key, m in (("3", ds["y"] == 0), ("4", ds["y"] == 1)):
            parts["Xs" + key].append(ds["F0"][m])
            parts["gs" + key].append(np.full(int(m.sum()), k))
    cat = {k: (np.vstack(v_) if k.startswith("X") else np.concatenate(v_)) for k, v_ in parts.items()}
    return {
        g: classifierTwoSample(cat["Xr" + key], cat["gr" + key], cat["Xs" + key], cat["gs" + key])
        for g, key in (("G3", "3"), ("G4", "4"))
    }


def withoutParticipants(bank: dict, drop: set, cfg: dict) -> dict:
    """The bank without some participants, with the artifact library rebuilt from the rest."""
    from .bank import buildLibrary

    keep = [e for e in bank["estimates"] if e["pid"] not in drop]
    return {
        "motion": [m for m in bank["motion"] if m["pid"] not in drop],
        "phys": [p for p in bank["phys"] if p["pid"] not in drop],
        "library": {**bank["library"], **buildLibrary(keep, cfg)},
        "source": bank["source"],
    }


def realismSplitHalf(real: list[Recording], bank: dict, cfg: dict, seed: int) -> dict:
    """Generator from one half of the participants, real windows of the other half, both ways."""
    pids = sorted(p["pid"] for p in bank["phys"])
    order = np.random.default_rng(seed).permutation(len(pids))
    halves = [{pids[i] for i in order[: len(pids) // 2]}, {pids[i] for i in order[len(pids) // 2 :]}]
    out = {"G3": {"halves": []}, "G4": {"halves": []}}
    for k, (gen, test) in enumerate(((halves[0], halves[1]), (halves[1], halves[0]))):
        sub = withoutParticipants(bank, test, cfg)
        r = realismPair([x for x in real if x.pid in test], sub, cfg, seed + k)
        for g in ("G3", "G4"):
            out[g]["halves"].append({"generatorPids": sorted(gen), "auc": r[g]["auc"]})
    for g in ("G3", "G4"):
        out[g]["auc"] = float(max(h["auc"] for h in out[g]["halves"]))
    return out


def realismPair(real: list[Recording], bank: dict, cfg: dict, seed: int) -> dict:
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
