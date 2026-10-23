"""Building synthetic data sets and evaluating them (docs/METHOD.md sections 5-8).

One *replicate* of one *condition* is a data set of synthetic participants whose
motion donors are the real participants (each used once) and whose physiology
donors are assigned by a derangement (nobody is paired with themself). Every
model and feature set is evaluated on it with leave-one-participant-out.
"""

from __future__ import annotations

import itertools
import json
import zlib
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

from .config import key
from .data import WIN
from .features import f0, fa, fa4, fb
from .generator import Condition, Generator, derangement
from .models import (
    auprc,
    auroc,
    bootstrapParticipants,
    losoPredict,
    sensitivityAtSpecificity,
    thresholdAtSpecificity,
)

featureSets = {
    "F0": ("F0",),
    "F0+FA": ("F0", "FA"),
    "F0+FA4": ("F0", "FA4"),
    "F0+FB": ("F0", "FB"),
}
deltas = {"dA": "F0+FA", "dA4": "F0+FA4", "dB": "F0+FB"}


# ------------------------------------------------------------------ bank
def loadBank(path: str | Path) -> dict:
    """Load the estimated parameters and motion donors written by ``cli estimate``."""
    path = Path(path)
    meta = json.loads(path.with_suffix(".json").read_text())
    if meta["library"].get("G") is None:  # written as null when no C4 gain could be estimated
        meta["library"]["G"] = float("nan")
    arrays = np.load(path.with_suffix(".npz"), allow_pickle=True)
    motion = []
    for m in meta["motion"]:
        pid = m["pid"]
        motion.append(
            {
                "pid": pid,
                "fs": m["fs"],
                "acc": arrays[f"acc{pid}"],
                "stage": arrays[f"stage{pid}"],
                "events": m["events"],
                "qAcc": m["qAcc"],
            }
        )
    return {"motion": motion, "phys": meta["phys"], "library": meta["library"], "source": meta["source"]}


# ------------------------------------------------------------------ data sets
def featuresFor(rec, cfg: dict) -> dict:
    fcfg = cfg["features"]
    fs = rec.fs
    return {
        "F0": f0(rec.y, fcfg["derivative"]),
        "FA": fa(rec.acc, fcfg["derivative"]),
        "FA4": fa4(rec.acc),
        "FB": fb(rec.y, rec.acc, fs, max(rec.meta["q"], 1e-6), max(rec.meta["qAcc"], 1e-6), **fcfg["fb"]),
    }


def makeDataset(
    bank: dict,
    cond: Condition,
    cfg: dict,
    seed: int,
    nParticipants: int | None = None,
    accShift: bool = False,
) -> dict:
    """Generate one replicate and its features.

    ``accShift`` circularly shifts each accelerometer record by a random lag of
    at least ``validation.shiftMinS`` *after* generation (negative control):
    movement stays real but no longer lines up with the artifacts.
    """
    rng = np.random.default_rng(seed)
    shiftRng = np.random.default_rng([seed, 1])
    gen = Generator(bank["library"], cfg)
    motion = bank["motion"]
    phys = bank["phys"]
    nAll = len(motion)
    order = np.arange(nAll)
    if nParticipants is not None and nParticipants < nAll:
        order = rng.choice(nAll, nParticipants, replace=False)
    perm = derangement(nAll, rng)
    recs, feats, ys, groups = [], [], [], []
    for gi, k in enumerate(order):
        m = motion[k]
        p = phys[int(perm[k])]
        rec = gen.generate(f"syn{gi:02d}", m, p, cond, rng)
        if accShift:  # separate stream, so the generated EDA and labels are identical to the unshifted data set
            lag = int(shiftRng.integers(int(cfg["validation"]["shiftMinS"] * m["fs"]), rec.acc.shape[0] - 1))
            rec.acc = np.roll(rec.acc, lag, axis=0)
        f = featuresFor(rec, cfg)
        nw = f["F0"].shape[0]
        recs.append(rec)
        feats.append(f)
        ys.append(rec.labels[:nw])
        groups.append(np.full(nw, gi))
    out = {name: np.vstack([f[name] for f in feats]) for name in ("F0", "FA", "FA4", "FB")}
    out.update(y=np.concatenate(ys), groups=np.concatenate(groups), records=recs)
    return out


def applyLabelNoise(y: np.ndarray, kappaTarget: float, rng: np.random.Generator) -> tuple[np.ndarray, float]:
    """Flip labels symmetrically until Cohen's kappa with the truth is about ``kappaTarget``.

    Simplest model of imperfect expert labelling (an assumption, tested separately).
    """
    lo, hi = 0.0, 0.5
    for iteration in range(30):
        f = 0.5 * (lo + hi)
        flip = rng.random(y.size) < f
        yn = np.where(flip, 1 - y, y)
        k = cohen_kappa_score(y, yn)
        if k > kappaTarget:
            lo = f
        else:
            hi = f
    return yn, float(k)


# ------------------------------------------------------------------ evaluation
def evaluate(
    ds: dict, cfg: dict, models: tuple = ("logreg",), seed: int = 0, labels: np.ndarray | None = None
) -> list[dict]:
    """Metrics for every model x feature set, plus the differences that form the map."""
    y = ds["y"] if labels is None else labels
    yTrue = ds["y"]
    g = ds["groups"]
    spec = cfg["analysis"]["specificity"]
    rows = []
    for model in models:
        preds = {}
        for name, parts in featureSets.items():
            X = np.hstack([ds[p] for p in parts])
            preds[name], chosen = losoPredict(X, y, g, model, cfg, seed)
        row = {"model": model, "nWindows": int(yTrue.size), "nPositive": int(yTrue.sum())}
        for name, s in preds.items():
            row[key("auprc", name)] = auprc(yTrue, s)
            row[key("auroc", name)] = auroc(yTrue, s)
            row[key("sens95", name)] = sensitivityAtSpecificity(yTrue, s, spec)
            lost = lostResponses(ds["records"], s, thresholdAtSpecificity(yTrue, s, spec))
            row[key("lostSpont", name)] = lost["spontaneous"]
            row[key("lostMotion", name)] = lost["motion"]
            for cls, frac in detectedByClass(ds["records"], s, thresholdAtSpecificity(yTrue, s, spec)).items():
                row[key("det", cls, name)] = frac
        for d, name in deltas.items():
            row[d] = row[key("auprc", name)] - row["auprcF0"]
            row[key(d, "auroc")] = row[key("auroc", name)] - row["aurocF0"]
        row["dBMinusDA4"] = row["dB"] - row["dA4"]
        row["predictions"] = preds
        rows.append(row)
    return rows


def recordOffsets(records):
    off, out = 0, []
    for r in records:
        out.append(off)
        off += r.labels.size
    return out


def lostResponses(records, scores: np.ndarray, thr: float) -> dict:
    """Share of true responses whose peak lies in a window flagged as artifact."""
    res = {"spontaneous": [], "motion": []}
    for r, off in zip(records, recordOffsets(records), strict=True):
        nw = r.labels.size
        for resp in r.responses:
            w = resp["peakSample"] // WIN
            if 0 <= w < nw and resp["amp"] >= 0.05:
                res[resp["origin"]].append(scores[off + w] > thr)
    return {k: (float(np.mean(v)) if v else float("nan")) for k, v in res.items()}


def detectedByClass(records, scores: np.ndarray, thr: float) -> dict:
    """Share of windows positive for each artifact class that the model flags."""
    hits: dict = {}
    for r, off in zip(records, recordOffsets(records), strict=True):
        nw = r.labels.size
        for cls, lab in r.meta["classWindows"].items():
            idx = np.where(lab[:nw] == 1)[0]
            hits.setdefault(cls, []).extend((scores[off + idx] > thr).tolist())
    return {k: float(np.mean(v)) for k, v in hits.items() if v}


# ------------------------------------------------------------------ grid
def gridConditions(cfg: dict) -> list[Condition]:
    gcfg = cfg["grid"]
    out = []
    for c, kappa, dmax, pir, comp in itertools.product(
        gcfg["c"], gcfg["kappa"], gcfg["deltaMax"], gcfg["piR"], gcfg["compositions"]
    ):
        out.append(
            Condition(c=c, kappa=kappa, deltaMax=dmax, piR=pir, classes=tuple(comp), prevalence=gcfg["prevalence"])
        )
    return out


def variantConditions(cfg: dict) -> list[Condition]:
    """Extra variants on part of the grid (c in {0, .5, 1}, kappa=1, deltaMax=1, piR=0)."""
    v = cfg["variants"]
    base = dict(kappa=1.0, deltaMax=1.0, piR=0.0)
    out = []
    for c in v["c"]:
        mix = tuple(cfg["grid"]["compositions"][-1])
        out.append(Condition(c=c, classes=("C2",), variant="V2", **base))
        out.append(Condition(c=c, classes=("C3",), variant="V2", **base))
        out.append(Condition(c=c, classes=mix, prevalence=v["highPrevalence"], **base))
        out.append(Condition(c=c, classes=mix, labelRule="majority", **base))
    out.append(Condition(c=1.0, classes=("C4",), **base))
    return out


def runCell(bank, cond: Condition, cfg: dict, rep: int, models: tuple, labelNoise: float | None = None) -> list[dict]:
    seed = cfg["seeds"]["grid"] * 1000003 + (zlib.crc32(cond.key().encode()) % 100000) * 101 + rep
    ds = makeDataset(bank, cond, cfg, seed)
    labels = None
    if labelNoise is not None:
        labels, unused = applyLabelNoise(ds["y"], labelNoise, np.random.default_rng(seed + 7))
    rows = evaluate(ds, cfg, models, seed, labels)
    for r in rows:
        r.pop("predictions")
        r.update(asdict(cond))
        r["classes"] = "+".join(cond.classes)
        r["rep"] = rep
        r["seed"] = seed
        r["prevalenceObserved"] = float(ds["y"].mean())
        r["nWindows"] = int(ds["y"].size)
        r["labelNoiseKappa"] = labelNoise
    return rows


def requiredReps(sdDelta: float, halfWidth: float, cap: int) -> int:
    """R = ceil((1.96 s / h)^2), limited to [2, cap] (docs/METHOD.md section 8)."""
    return int(min(cap, max(2, np.ceil((1.96 * sdDelta / halfWidth) ** 2))))


def summarise(df: pd.DataFrame, by: list[str], value: str, rule: dict) -> pd.DataFrame:
    """Mean, 95 % t-interval and the decision rule per cell."""
    from scipy import stats

    def agg(x):
        x = x.dropna().to_numpy()
        n = x.size
        m = x.mean() if n else np.nan
        h = stats.t.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan
        return pd.Series(
            {
                "mean": m,
                "lo": m - h,
                "hi": m + h,
                "n": n,
                "useful": bool(m >= rule["minDelta"] and (m - h) > 0) if n > 1 else False,
            }
        )

    return df.groupby(by)[value].apply(agg).unstack()


def compatibility(df: pd.DataFrame, observed: float, by: list[str], value: str = "dA4") -> pd.DataFrame:
    """Goal 5: is the observed real difference inside the central 95 % of a cell's replicates?"""

    def f(x):
        lo, hi = np.percentile(x, [2.5, 97.5])
        return pd.Series({"p2.5": lo, "p97.5": hi, "compatible": bool(lo <= observed <= hi), "n": x.size})

    return df.groupby(by)[value].apply(f).unstack()


def bootstrapPower(bank, cond: Condition, cfg: dict, nPart: int, reps: int, models=("logreg",)) -> dict:
    """Goal 6: share of replicates whose participant-bootstrap CI of dA4 and dB excludes 0."""
    res = {"dA4": [], "dB": []}
    for rep in range(reps):
        seed = cfg["seeds"]["power"] * 7919 + nPart * 1009 + rep
        ds = makeDataset(bank, cond, cfg, seed, nParticipants=nPart)
        row = evaluate(ds, cfg, models, seed)[0]
        p = row["predictions"]
        for d, name in (("dA4", "F0+FA4"), ("dB", "F0+FB")):
            lo, hi, unused = bootstrapParticipants(
                ds["y"], p["F0"], p[name], ds["groups"], cfg["analysis"]["bootstrapN"], seed
            )
            res[d].append(lo > 0)
    return {d: float(np.mean(v)) for d, v in res.items()}
