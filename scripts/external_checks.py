"""Checks on independent data, as planned in docs/EXTERNAL-CHECKS.md.

Usage:
    python scripts/external_checks.py edabe <edabe-8hz.npz>  [--reps 10]
    python scripts/external_checks.py wesad <wesad-wrist.npz> [--reps 50] [--jobs 2]

Writes results/external-edabe.json or results/external-wesad.json (and a CSV of the
WESAD replicates). The bank of `estimate` (results/bank-utd) and, for the UTD-trained
detector, the raw UTD recordings are needed.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from edasynth import cli
from edasynth.config import loadConfig
from edasynth.data import WIN, Recording
from edasynth.estimation import motionEvents, quantisationStep
from edasynth.experiment import evaluate, loadBank, makeDataset
from edasynth.features import f0
from edasynth.generator import Condition
from edasynth.models import auprc, losoPredict, makeModel
from edasynth.validation import classifierTwoSample, tCi

MIX = ("C1", "C2", "C3", "C5")


def setup():
    cfg = loadConfig(None)
    cfg["repoRoot"] = str(Path(cli.__file__).resolve().parents[2])
    res = cli.resultsDir(cfg)
    return cfg, res, loadBank(res / "bank-utd", cfg)


def zWithin(X: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Each feature z-scored within each recording (constant features become 0)."""
    Z = np.zeros_like(X, dtype=float)
    for k in np.unique(g):
        m = g == k
        sd = X[m].std(0)
        Z[m] = (X[m] - X[m].mean(0)) / np.where(sd > 0, sd, 1.0)
    return Z


# ------------------------------------------------------------------ EDABE
def loadEdabe(path: str, deriv: str):
    """F0 features, window labels (any labelled sample) and recording index of every EDABE window."""
    z = np.load(path)
    names = sorted({k.rsplit("__", 1)[0] for k in z.files})
    X, y, g, info = [], [], [], []
    for k, name in enumerate(names):
        raw = np.asarray(z[name + "__raw"], float)
        lab = np.asarray(z[name + "__label"], int)
        nz = np.flatnonzero(raw != 0)  # leading/trailing zeros: device not yet / no longer recording
        if nz.size == 0:
            continue
        raw, lab = raw[nz[0] : nz[-1] + 1], lab[nz[0] : nz[-1] + 1]
        nw = min(raw.size, lab.size) // WIN
        F = f0(raw[: nw * WIN], deriv)[:nw]
        yy = lab[: nw * WIN].reshape(nw, WIN).max(1)
        ok = np.isfinite(F).all(1)
        X.append(F[ok])
        y.append(yy[ok])
        g.append(np.full(int(ok.sum()), k))
        info.append({"name": name, "windows": int(nw), "dropped": int((~ok).sum()), "artifactShare": float(yy.mean())})
    return np.vstack(X), np.concatenate(y), np.concatenate(g), info


def synthetic(bank, cfg, seed):
    cond = Condition(c=0.5, kappa=1.0, deltaMax=1.0, piR=0.0, classes=MIX, prevalence=cfg["grid"]["prevalence"])
    ds = makeDataset(bank, cond, cfg, seed)
    return ds["F0"], ds["y"], ds["groups"]


def utdTraining(cfg):
    deriv = cfg["features"]["derivative"]
    X, y, g = [], [], []
    for k, rec in enumerate(cli.utdAligned(cfg)):
        F = f0(rec.eda, deriv)
        nw = min(F.shape[0], rec.labels3.shape[0])
        X.append(F[:nw])
        y.append(rec.majority()[:nw])
        g.append(np.full(nw, k))
    return np.vstack(X), np.concatenate(y), np.concatenate(g)


def fitScore(Xtr, ytr, Xte, cfg):
    m = makeModel("logreg", cfg, C=1.0).fit(Xtr, ytr)
    return m.predict_proba(Xte)[:, 1]


def cmdEdabe(args):
    cfg, res, bank = setup()
    deriv = cfg["features"]["derivative"]
    Xe, ye, ge, info = loadEdabe(args.data, deriv)
    Ze = zWithin(Xe, ge)
    out = {
        "recordings": len(info),
        "windows": int(ye.size),
        "artifactWindowShare": float(ye.mean()),
        "perRecording": info,
        "reps": args.reps,
    }
    rows = []
    for r in range(args.reps):
        Xs, ys, gs = synthetic(bank, cfg, 777000 + r)
        Zs = zWithin(Xs, gs)
        row = {"rep": r}
        for tag, A, B in (("z", Ze, Zs), ("raw", Xe, Xs)):
            row[f"E1{tag}"] = classifierTwoSample(A[ye == 1], ge[ye == 1], B[ys == 1], gs[ys == 1])["auc"]
            row[f"E2{tag}"] = classifierTwoSample(A[ye == 0], ge[ye == 0], B[ys == 0], gs[ys == 0])["auc"]
        row["E3a"] = auprc(ye, fitScore(Zs, ys, Ze, cfg))
        rows.append(row)
        cli.log(f"EDABE rep {r}: {row}")
    df = pd.DataFrame(rows)
    Xu, yu, gu = utdTraining(cfg)
    e3b = auprc(ye, fitScore(zWithin(Xu, gu), yu, Ze, cfg))
    e3c = auprc(ye, losoPredict(Ze, ye, ge, "logreg", cfg)[0])
    lim = cfg["validation"]["c2stAucMax"]
    for k in ("E1z", "E1raw", "E2z", "E2raw", "E3a"):
        out[k] = {"median": float(df[k].median()), "min": float(df[k].min()), "max": float(df[k].max())}
    out["E1z"]["passed"] = bool(out["E1z"]["median"] <= lim)
    out["E2z"]["passed"] = bool(out["E2z"]["median"] <= lim)
    out["E3b"] = float(e3b)
    out["E3c"] = float(e3c)
    out["E3passed"] = bool(out["E3a"]["median"] >= e3b - 0.05)
    out["replicates"] = rows
    (res / "external-edabe.json").write_text(json.dumps(out, indent=1))
    cli.log(json.dumps({k: v for k, v in out.items() if k not in ("perRecording", "replicates")}, indent=1))


# ------------------------------------------------------------------ WESAD
def wesadMotion(path: str, cfg: dict) -> list[dict]:
    """WESAD wrist recordings as movement donors (docs/EXTERNAL-CHECKS.md, W1)."""
    z = np.load(path)
    subjects = sorted({k.split("__")[0] for k in z.files}, key=lambda s: int(s.lstrip("S")))
    motion = []
    for s in subjects:
        acc32 = np.asarray(z[s + "__acc"], float) / 64.0
        lab4 = np.asarray(z[s + "__label4hz"], int)
        n8 = acc32.shape[0] // 4
        acc = acc32[: n8 * 4].reshape(n8, 4, 3).mean(1)
        lab8 = np.repeat(lab4, 2)[:n8]
        if lab8.size < n8:
            lab8 = np.r_[lab8, np.zeros(n8 - lab8.size, int)]
        idx = np.where((lab8 >= 1) & (lab8 <= 4))[0]
        a, b = int(idx[0]), int(idx[-1]) + 1
        acc, lab8 = acc[a:b], lab8[a:b]
        stage = np.where(np.isin(lab8, [1, 4]), "rest", "load")
        rec = Recording(pid=f"W{s}", fs=8.0, eda=np.zeros(acc.shape[0]), acc=acc, stage=stage)
        events, theta = motionEvents(rec, cfg)
        motion.append(
            {
                "pid": rec.pid,
                "fs": 8.0,
                "acc": acc,
                "stage": stage,
                "events": events,
                "qAcc": quantisationStep(acc[:, 0]),
                "durationMin": acc.shape[0] / 8.0 / 60.0,
                "eventsPerMin": len(events) / (acc.shape[0] / 8.0 / 60.0),
            }
        )
    return motion


def wesadReplicate(bank, cond, cfg, seed):
    ds = makeDataset(bank, cond, cfg, seed)
    row = evaluate(ds, cfg, ("logreg",), seed)[0]
    row.pop("predictions", None)
    return {k: row[k] for k in ("dA", "dA4", "dB", "auprcF0", "nWindows", "nPositive")}


def cmdWesad(args):
    from joblib import Parallel, delayed

    cfg, res, bank = setup()
    motion = wesadMotion(args.data, cfg)
    utdRate = [len(m["events"]) / (m["acc"].shape[0] / m["fs"] / 60.0) for m in bank["motion"]]
    out = {
        "donors": [{k: m[k] for k in ("pid", "durationMin", "eventsPerMin", "qAcc")} for m in motion],
        "utdEventsPerMinMedian": float(np.median(utdRate)),
        "wesadEventsPerMinMedian": float(np.median([m["eventsPerMin"] for m in motion])),
    }
    wbank = {**bank, "motion": [{k: m[k] for k in ("pid", "fs", "acc", "stage", "events", "qAcc")} for m in motion]}
    grid = cli.readGrid(res / "grid-metrics-utd.csv", cfg)
    rows, cells = [], []
    part = res / "external-wesad.json"
    if part.exists():  # resume: cells already finished by an interrupted run are kept
        prev = json.loads(part.read_text())
        cells = prev.get("cells", [])
        rows = pd.read_csv(res / "external-wesad-replicates.csv").to_dict("records")
    for c in (0.0, 0.5, 1.0):
        if any(x["c"] == c for x in cells):
            continue
        cond = Condition(c=c, kappa=1.0, deltaMax=1.0, piR=0.0, classes=MIX, prevalence=cfg["grid"]["prevalence"])
        reps = Parallel(n_jobs=args.jobs)(
            delayed(wesadReplicate)(wbank, cond, cfg, 888000 + int(c * 100) * 1000 + r) for r in range(args.reps)
        )
        for r, row in enumerate(reps):
            rows.append({"c": c, "rep": r, **row})
        u = grid[
            (grid.classes == "+".join(MIX))
            & (grid.c == c)
            & (grid.kappa == 1.0)
            & (grid.deltaMax == 1.0)
            & (grid.piR == 0.0)
        ]
        cell = {"c": c}
        for d in ("dA4", "dB"):
            for tag, x in (("wesad", [row[d] for row in reps]), ("utd", u[d].to_numpy())):
                m, lo, hi = tCi(x)
                cell[f"{d}_{tag}"] = {
                    "mean": m,
                    "lo": lo,
                    "hi": hi,
                    "n": int(len(x)),
                    "useful": bool(m >= 0.01 and lo > 0),
                }
        cell["sameDecisionFB"] = cell["dB_wesad"]["useful"] == cell["dB_utd"]["useful"]
        cells.append(cell)
        cli.log(f"WESAD c={c}: {cell}")
        pd.DataFrame(rows).to_csv(res / "external-wesad-replicates.csv", index=False)  # kept if a later cell is cut
        (res / "external-wesad.json").write_text(json.dumps({**out, "cells": cells, "partial": True}, indent=1))
    out["cells"] = cells
    out["passed"] = bool(all(c["sameDecisionFB"] for c in cells))
    pd.DataFrame(rows).to_csv(res / "external-wesad-replicates.csv", index=False)
    (res / "external-wesad.json").write_text(json.dumps(out, indent=1))
    cli.log(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["edabe", "wesad"])
    ap.add_argument("data")
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    if args.which == "edabe":
        args.reps = args.reps or 10
        cmdEdabe(args)
    else:
        args.reps = args.reps or 50
        cmdWesad(args)


if __name__ == "__main__":
    main()
