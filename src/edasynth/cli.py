"""Command-line entry point: ``python -m edasynth <command>`` (see ``--help``).

Commands follow the order of the method; ``make all`` runs them in sequence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .config import key, loadConfig, resolve

resultsDefault = "results"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def writeJson(path: Path, obj):
    from .bank import toJson

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(toJson(obj))


def resultsDir(cfg) -> Path:
    p = resolve(cfg, "results")
    p.mkdir(parents=True, exist_ok=True)
    return p


# ------------------------------------------------------------------ data
def cmdCheckData(cfg, args):
    """Verify SHA-256 checksums of the downloaded files."""
    bad = 0
    lists = (
        ("noneeg", resolve(cfg, "noneeg") / "SHA256SUMS.txt"),  # official list shipped with the data
        ("zhang", Path(cfg["repoRoot"]) / "data/checksums/zhang-eda-artifact-detection-7c4bd46.sha256"),
    )
    for source, listfile in lists:
        root = resolve(cfg, source)
        if not listfile.exists():
            print(f"MISSING  checksum list {listfile}")
            bad += 1
            continue
        lines = listfile.read_text().splitlines()
        for ln in lines:
            if not ln.strip():
                continue
            digest, name = ln.split(maxsplit=1)
            f = root / name
            if not f.exists():
                print(f"MISSING  {source}/{name}")
                bad += 1
                continue
            got = hashlib.sha256(f.read_bytes()).hexdigest()
            ok = got == digest
            bad += not ok
            if not ok or args.verbose:
                print(f"{'OK      ' if ok else 'MISMATCH'} {source}/{name}")
    print("all files verified" if bad == 0 else f"{bad} problem(s)")
    return 1 if bad else 0


def utd(cfg, withLabels=True):
    from .data import loadUtdLabels, loadUtdParticipant

    labels = loadUtdLabels(resolve(cfg, "zhang")) if withLabels else None
    excl = set(cfg["data"]["excludeSubjects"])
    recs = []
    for s in range(1, 21):
        if s in excl:
            continue
        recs.append(loadUtdParticipant(resolve(cfg, "noneeg"), s, labels))
    return recs


def utdAligned(cfg):
    """UTD recordings accepted by `align`, with the alignment offset applied."""
    path = resultsDir(cfg) / "alignment.csv"
    if not path.exists():
        raise SystemExit("results/alignment.csv not found: run `python -m edasynth align` first")
    align = pd.read_csv(path)
    ok = align[align.accepted]
    off = dict(zip(ok.pid, ok.offsetSamples, strict=True))
    return [r.cropped(int(off[r.pid])) for r in utd(cfg) if r.pid in off]


def cmdInspect(cfg, args):
    """Units, ranges, annotation codes and stage parsing for every UTD record."""
    from .data import stageType

    rows = []
    for r in utd(cfg):
        mag = np.sqrt((r.acc**2).sum(1))
        rest = np.array([stageType(str(s)) == "rest" for s in r.stage]) if r.stage is not None else np.ones(r.n, bool)
        rows.append(
            {
                "pid": r.pid,
                "nSamples": r.n,
                "nWindows": r.nWindows,
                "fs": r.fs,
                "units": r.units,
                "edaMin": float(r.eda.min()),
                "edaMedian": float(np.median(r.eda)),
                "edaMax": float(r.eda.max()),
                "accMagMedianRest": float(np.median(mag[rest])),
                "nAnnotations": len(r.meta.get("annotations", [])),
                "annotationCodes": sorted({a[1] for a in r.meta.get("annotations", [])}),
                "annotationAux": sorted({a[2] for a in r.meta.get("annotations", [])}),
                "stageError": r.meta.get("stageError", ""),
            }
        )
    out = resultsDir(cfg) / "inspect-utd.json"
    writeJson(out, rows)
    df = pd.DataFrame(rows)
    print(
        df[
            [
                "pid",
                "nSamples",
                "edaMin",
                "edaMedian",
                "edaMax",
                "accMagMedianRest",
                "nAnnotations",
                "stageError",
            ]
        ].to_string()
    )
    log(f"written {out}")


def cmdAlign(cfg, args):
    """Check that the expert labels line up with the raw signals (docs/METHOD.md section 2.4)."""
    from .alignment import checkParticipant
    from .data import loadUtdLabels, loadUtdStoredFeatures

    labels = loadUtdLabels(resolve(cfg, "zhang"))
    edaStored, unused = loadUtdStoredFeatures(resolve(cfg, "zhang"))
    rows = []
    for r in utd(cfg, withLabels=False):
        g = int(r.pid[1:]) - 1
        m = labels.group.to_numpy() == g
        res = checkParticipant(r.eda, edaStored[m][:, :4], cfg["data"]["alignmentMinR"])
        res.update(pid=r.pid, group=g, artifactShare=float(labels.majority[m].mean()))
        rows.append(res)
        log(
            f"{r.pid}: windows raw/labels {res['nWindowsRaw']}/{res['nWindowsLabels']}, "
            f"rMin={res['rMin']:.4f}, shift={res['shiftWindows']}w+{res['shiftSamples']}s, "
            f"{'ACCEPTED' if res['accepted'] else 'REJECTED'}"
        )
    df = pd.DataFrame(rows)
    out = resultsDir(cfg) / "alignment.csv"
    df.to_csv(out, index=False)
    if args.plots:
        import matplotlib.pyplot as plt

        from .figures import CB, saveFigure

        offsets = dict(zip(df.pid, df.offsetSamples, strict=True))
        for r0 in utd(cfg):
            r = r0.cropped(int(offsets.get(r0.pid, 0)))
            fig, ax = plt.subplots(3, 1, figsize=(12, 5), sharex=True)
            t = np.arange(r.n) / r.fs / 60
            ax[0].plot(t, r.eda, lw=0.6, color=CB[0])
            ax[0].set_ylabel("EDA")
            ax[1].plot(t, np.sqrt((r.acc**2).sum(1)), lw=0.5, color=CB[1])
            ax[1].set_ylabel("|a|, g")
            for k in range(3):
                ax[2].fill_between(
                    t[: r.labels3.shape[0] * 40], k, k + np.repeat(r.labels3[:, k], 40) * 0.8, step="pre", color=CB[3]
                )
            ax[2].set_ylabel("experts 1-3")
            ax[2].set_xlabel("time, min")
            if r.stage is not None:
                ch = np.r_[0, np.where(r.stage[1:] != r.stage[:-1])[0] + 1]
                for c in ch:
                    for a in ax:
                        a.axvline(t[c], color="0.5", lw=0.6, ls=":")
            fig.suptitle(r.pid)
            saveFigure(fig, resolve(cfg, "figures") / "alignment", f"step2-{r.pid}")
    log(f"written {out}; accepted {int(df.accepted.sum())}/{len(df)}")


# ------------------------------------------------------------------ real data
def cmdRealdata(cfg, args):
    """Baseline on real UTD: stored features (reproduction) and/or raw recomputation (G7)."""
    from .data import loadUtdLabels, loadUtdStoredFeatures
    from .features import f0, fa4
    from .realdata import compare, shuffledControl, singleFeatureAuprc

    res = {}
    labels = loadUtdLabels(resolve(cfg, "zhang"))
    y = labels.majority.to_numpy()
    g = labels.group.to_numpy()
    rcfg = realdataCfg(cfg)
    if args.stored:
        eda, acc = loadUtdStoredFeatures(resolve(cfg, "zhang"))
        log("stored features: logistic regression, LOSO")
        res["storedLogreg"] = compare(eda, np.hstack([eda, acc[:, :4]]), y, g, rcfg, "logreg")
        log("stored features: random forest, LOSO")
        res["storedRf"] = compare(eda, np.hstack([eda, acc[:, :4]]), y, g, rcfg, "rf")
        log("stored features: shuffled-acceleration control")
        res["storedShuffled"] = shuffledControl(eda, acc[:, :4], y, g, rcfg, cfg["realdata"]["nShuffles"], 0)
        res["storedSingleFeatureAuprc"] = singleFeatureAuprc(eda, y)
    if args.raw:
        from scipy import stats

        accRows, edaRows, yy, gg = [], [], [], []
        for r in utdAligned(cfg):
            nw = r.labels3.shape[0]
            # z-scored within participant, as in the published feature matrices
            edaRows.append(np.nan_to_num(stats.zscore(f0(r.eda, cfg["features"]["derivative"])[:nw], axis=0)))
            accRows.append(np.nan_to_num(stats.zscore(fa4(r.acc)[:nw], axis=0)))
            yy.append(r.majority())
            gg.append(np.full(nw, int(r.pid[1:])))
        E, A = np.vstack(edaRows), np.vstack(accRows)
        Y, G = np.concatenate(yy), np.concatenate(gg)
        res["rawLogreg"] = compare(E, np.hstack([E, A]), Y, G, rcfg, "logreg")
        from .validation import g7RealConsistency

        res["G7"] = g7RealConsistency(res["rawLogreg"], cfg)
    for v in res.values():
        if isinstance(v, dict):
            v.pop("scores0", None)
            v.pop("scores1", None)
    out = resultsDir(cfg) / ("realdata.json" if not args.out else args.out)
    writeJson(out, res)
    for k, v in res.items():
        if isinstance(v, dict) and "delta" in v:
            log(
                f"{k}: AUPRC {v['auprc0']:.4f} -> {v['auprc1']:.4f}, delta {v['delta']:+.4f}, "
                f"95% CI [{v['ci95'][0]:+.4f}, {v['ci95'][1]:+.4f}], Wilcoxon p={v['wilcoxonP']:.3f}"
            )
    log(f"written {out}")


def realdataCfg(cfg):
    from .config import merge

    return merge(cfg, {"models": cfg["realdata"]["models"]})


# ------------------------------------------------------------------ estimation
def cmdEstimate(cfg, args):
    """Estimate generator parameters (docs/METHOD.md section 3) and write the bank."""
    from .bank import buildBank

    if args.source == "utd":
        recs = utdAligned(cfg)
        bad = [r.pid for r in recs if r.stage is None]
        if bad:
            raise SystemExit(f"stages could not be parsed for {bad}; run `inspect` first")
    else:
        from .data import loadAww

        recs = loadAww(resolve(cfg, "zhang"))
        log("AWW is ONE person: donors are different recording days, not different people (development only)")
    out = resultsDir(cfg) / f"bank-{args.source}"
    meta = buildBank(recs, cfg, args.source, out, log, nJobs=args.jobs)
    lib = meta["library"]
    log(f"episodes by class: {lib['episodeCounts']}; C4 gain G = {lib['G']:.4f}")
    for e in meta.get("estimates", []):
        for w in e.get("warnings", []):
            log(f"WARNING {e['pid']}: {w}")
    log(f"written {out}.json / .npz")


def loadSourceBank(cfg, args):
    from .experiment import loadBank

    return loadBank(resultsDir(cfg) / f"bank-{args.source}", cfg)


# ------------------------------------------------------------------ detectors (goal 2)
def cmdDetectors(cfg, args):
    from .detectors import runAll
    from .experiment import makeDataset
    from .generator import Condition
    from .models import auprc

    bank = loadSourceBank(cfg, args)
    rows = []
    for ci, cls in enumerate(cfg["detectors"]["classes"]):
        for rep in range(args.reps):
            cond = Condition(c=0.5, classes=(cls,), prevalence=cfg["grid"]["prevalence"])
            ds = makeDataset(bank, cond, cfg, cfg["seeds"]["detectors"] + rep * 31 + ci * 1000)
            flags = {}
            for r in ds["records"]:
                for name, (kind, val) in runAll(r.y, r.fs, cfg).items():
                    flags.setdefault(name, (kind, []))[1].append(val[: r.labels.size])
            y = ds["y"]
            if y.sum() == 0:
                log(f"{cls}: no positive windows generated (e.g. C4 gain ~0); sensitivity not defined")
                rows += [
                    {
                        "class": cls,
                        "rep": rep,
                        "detector": name,
                        "sensitivity": float("nan"),
                        "specificity": float("nan"),
                        "auprc": float("nan"),
                        "nPositive": 0,
                    }
                    for name in flags
                ]
                continue
            for name, (kind, parts) in flags.items():
                s = np.concatenate(parts).astype(float)
                if kind == "flag":
                    rows.append(
                        {
                            "class": cls,
                            "rep": rep,
                            "detector": name,
                            "sensitivity": float(s[y == 1].mean()),
                            "specificity": float(1 - s[y == 0].mean()),
                            "auprc": float("nan"),
                            "nPositive": int(y.sum()),
                        }
                    )
                else:
                    from .models import sensitivityAtSpecificity

                    rows.append(
                        {
                            "class": cls,
                            "rep": rep,
                            "detector": name,
                            "sensitivity": sensitivityAtSpecificity(y, s, 0.95),
                            "specificity": 0.95,
                            "auprc": auprc(y, s),
                            "nPositive": int(y.sum()),
                        }
                    )
        log(f"detectors done for {cls}")
    df = pd.DataFrame(rows)
    out = resultsDir(cfg) / "detectors-by-class.csv"
    df.to_csv(out, index=False)
    print(df.groupby(["detector", "class"]).sensitivity.mean().unstack().round(3).to_string())
    log(f"written {out}")


# ------------------------------------------------------------------ validation
def cmdValidate(cfg, args):
    from . import validation as V

    bank = loadSourceBank(cfg, args)
    rep = {"source": args.source, "version": __version__}
    if args.pytest:
        r = subprocess.run([sys.executable, "-m", "pytest", "-q"], capture_output=True, text=True)
        rep["G1"] = {"passed": r.returncode == 0, "tail": r.stdout.strip().splitlines()[-1:]}
    seed = cfg["seeds"]["validation"]
    out = resultsDir(cfg) / f"validation-{args.source}.json"

    def save():  # written after every gate, so a crash later does not lose finished gates
        writeJson(out, rep)

    if "G2" in args.gates:
        log("G2 parameter recovery")
        rep["G2"] = V.g2ParameterRecovery(bank, cfg, cfg["validation"]["g2Participants"], seed)
        save()
    if "G3" in args.gates or "G4" in args.gates:
        if args.source == "utd":
            log("G3/G4 realism")
            rep.update(V.g3G4Realism(utdAligned(cfg), bank, cfg, seed))
        else:
            from .data import loadAww

            log("G3/G4 realism against AWW (development only)")
            rep.update(V.g3G4Realism(loadAww(resolve(cfg, "zhang")), bank, cfg, seed))
        save()
    if "G5" in args.gates:
        log("G5 negative controls")
        rep["G5"] = V.g5NegativeControls(bank, cfg, seed + 1)
        save()
    if "G6" in args.gates:
        log("G6 positive controls")
        rep["G6"] = V.g6PositiveControls(bank, cfg, seed + 2)
        save()
    rj = resultsDir(cfg) / "realdata.json"
    if "G7" in args.gates and rj.exists():
        rd = json.loads(rj.read_text())
        if "G7" in rd:
            rep["G7"] = rd["G7"]
    save()
    for k, v in rep.items():
        if isinstance(v, dict) and "passed" in v:
            log(f"{k}: {'PASS' if v['passed'] else 'FAIL'}")
    log(f"written {out}")


# ------------------------------------------------------------------ experiment
def runJobs(cfg, bank, jobs, nJobs, outCsv, models, labelNoise=None):
    from joblib import Parallel, delayed

    from .experiment import runCell

    done = set()
    if outCsv.exists():
        old = pd.read_csv(outCsv)
        ids = set(old["bankId"].astype(str)) if "bankId" in old else {"none"}
        if ids != {bank["bankId"]}:
            raise SystemExit(
                f"{outCsv} holds results made with another parameter bank ({', '.join(sorted(ids))}; "
                f"current {bank['bankId']}). Move or delete it and run the command again."
            )
        done = set(zip(old["key"], old["rep"], strict=True))
    todo = [(c, r) for c, r in jobs if (c.key(), r) not in done]
    log(f"{len(todo)} replicate-cells to run ({len(jobs) - len(todo)} already done)")
    batch = max(1, nJobs * 2)
    for i in range(0, len(todo), batch):
        chunk = todo[i : i + batch]
        res = Parallel(n_jobs=nJobs)(delayed(runCell)(bank, c, cfg, r, models, labelNoise) for c, r in chunk)
        rows = []
        for (c, rep), rr in zip(chunk, res, strict=True):
            for row in rr:
                row["key"] = c.key()
                row["bankId"] = bank["bankId"]
                row["version"] = __version__
                rows.append(row)
        df = pd.DataFrame(rows)
        # rows of different cells can have different columns (e.g. detC2* only for
        # mixtures), so the file is rewritten with aligned columns instead of appended
        if outCsv.exists():
            df = pd.concat([pd.read_csv(outCsv), df], ignore_index=True, sort=False)
        tmp = outCsv.with_suffix(".tmp")
        df.to_csv(tmp, index=False)
        tmp.replace(outCsv)
        log(f"{min(i + batch, len(todo))}/{len(todo)}")


def cmdPilot(cfg, args):
    from .experiment import requiredReps
    from .generator import Condition

    bank = loadSourceBank(cfg, args)
    mix = tuple(cfg["grid"]["compositions"][-1])
    conds = [Condition(c=c, kappa=k, classes=mix) for c in (0.0, 1.0) for k in (0.0, 1.0)]
    jobs = [(c, r) for c in conds for r in range(cfg["pilot"]["reps"])]
    out = resultsDir(cfg) / f"pilot-{args.source}.csv"
    runJobs(cfg, bank, jobs, args.jobs, out, ("logreg",))
    df = pd.read_csv(out)
    sd = float(df.groupby("key")[["dA4", "dB"]].std().max().max())
    R = requiredReps(sd, cfg["pilot"]["halfWidth"], cfg["pilot"]["maxReps"])
    writeJson(resultsDir(cfg) / f"pilot-{args.source}-R.json", {"maxSdDelta": sd, "R": R})
    log(f"largest SD of delta across pilot cells = {sd:.4f} -> R = {R}")


def cmdGrid(cfg, args):
    from .experiment import gridToRun, variantConditions

    bank = loadSourceBank(cfg, args)
    R = args.reps
    if R is None:
        R = json.loads((resultsDir(cfg) / f"pilot-{args.source}-R.json").read_text())["R"]
    conds = variantConditions(cfg) if args.variants else gridToRun(cfg)
    models = tuple(args.models.split(","))
    jobs = [(c, r) for c in conds for r in range(R)]
    name = "variants" if args.variants else "grid"
    if args.labelNoise:
        name += "-labelnoise"
    out = resultsDir(cfg) / f"{name}-metrics-{args.source}.csv"
    runJobs(cfg, bank, jobs, args.jobs, out, models, cfg["analysis"]["labelNoiseKappa"] if args.labelNoise else None)
    log(f"written {out}")


def cmdAnalyze(cfg, args):
    from .experiment import compatibility, summarise

    res = resultsDir(cfg)
    df = readGrid(res / f"grid-metrics-{args.source}.csv", cfg)
    by = ["c", "kappa", "deltaMax", "piR", "classes"]
    parts = []
    for d in ("dA", "dA4", "dB", "dBMinusDA4"):
        s = summarise(df, by, d, cfg["analysis"]["decision"])
        s.columns = [key(c, d) for c in s.columns]
        parts.append(s)
    summ = pd.concat(parts, axis=1)
    summ.to_csv(res / f"grid-summary-{args.source}.csv")
    comp = compatibility(
        df[df.classes == "+".join(cfg["grid"]["compositions"][-1])],
        cfg["analysis"]["observedDelta"],
        ["c", "kappa", "deltaMax", "piR"],
    )
    comp.to_csv(res / f"compatibility-{args.source}.csv")
    ln = res / f"grid-labelnoise-metrics-{args.source}.csv"
    if ln.exists():  # goal 5 with noisy training labels
        dn = readGrid(ln, cfg)
        compN = compatibility(
            dn[dn.classes == "+".join(cfg["grid"]["compositions"][-1])],
            cfg["analysis"]["observedDelta"],
            ["c", "kappa", "deltaMax", "piR"],
        )
        compN.to_csv(res / f"compatibility-labelnoise-{args.source}.csv")
        log(f"compatible cells with label noise: {int(compN.compatible.sum())}/{len(compN)}")
    ref = cfg["analysis"]["sensitivityReference"]
    mix = "+".join(cfg["grid"]["compositions"][-1])
    sens = []
    base = df[df.classes == mix]
    for p, levels in (
        ("c", cfg["grid"]["c"]),
        ("kappa", cfg["grid"]["kappa"]),
        ("deltaMax", cfg["grid"]["deltaMax"]),
        ("piR", cfg["grid"]["piR"]),
    ):
        sub = base.copy()
        for q, val in ref.items():
            if q != p:
                sub = sub[np.isclose(sub[q], val)]
        means = sub.groupby(p)[["dA4", "dB"]].mean()
        sens.append(
            {
                "parameter": p,
                "levels": len(levels),
                "rangeDA4": float(means.dA4.max() - means.dA4.min()),
                "rangeDB": float(means.dB.max() - means.dB.min()),
            }
        )
    pd.DataFrame(sens).sort_values("rangeDB", ascending=False).to_csv(
        res / f"sensitivity-{args.source}.csv", index=False
    )
    log(f"compatible cells: {int(comp.compatible.sum())}/{len(comp)}")
    print(pd.DataFrame(sens).to_string())


def cmdPower(cfg, args):
    from .experiment import bootstrapPower
    from .generator import Condition

    bank = loadSourceBank(cfg, args)
    res = resultsDir(cfg)
    comp = pd.read_csv(res / f"compatibility-{args.source}.csv")
    mix = tuple(cfg["grid"]["compositions"][-1])
    gm = readGrid(res / f"grid-metrics-{args.source}.csv", cfg)
    gm = gm[gm.classes == "+".join(mix)]
    by = ["c", "kappa", "deltaMax", "piR"]
    comp = comp[comp.compatible].merge(gm.groupby(by, as_index=False)["dB"].mean(), on=by)
    comp = selectPowerConditions(comp, cfg["power"]["nConditions"])
    rows = []
    for rowIndex, c in comp.iterrows():
        cond = Condition(c=c.c, kappa=c.kappa, deltaMax=c.deltaMax, piR=c.piR, classes=mix)
        for n in cfg["power"]["participants"]:
            pw = bootstrapPower(bank, cond, cfg, n, cfg["power"]["reps"])
            rows.append({"condition": cond.key(), "nParticipants": n, **{key("power", k): v for k, v in pw.items()}})
            log(f"{cond.key()} n={n}: {pw}")
    pd.DataFrame(rows).to_csv(res / f"protocol-planning-{args.source}.csv", index=False)


def readGrid(path: Path, cfg: dict) -> pd.DataFrame:
    """Logistic-regression rows of a grid file, with the equivalent cells filled in."""
    from .experiment import expandGrid

    df = pd.read_csv(path)
    return expandGrid(df[df.model == "logreg"], cfg)


def selectPowerConditions(comp: pd.DataFrame, k: int) -> pd.DataFrame:
    """k compatible cells spread evenly over their mean gain dB (docs/METHOD.md section 8).

    Sorted by mean dB, the cells at the quantiles 1/(k+1), ..., k/(k+1) are taken, so the
    power is computed for small, middle and large gains that agree with the real data
    instead of the first cells in grid order [decision].
    """
    comp = comp.sort_values("dB").reset_index(drop=True)
    if len(comp) <= k:
        return comp
    idx = sorted({int(round(q * (len(comp) - 1))) for q in np.arange(1, k + 1) / (k + 1)})
    return comp.iloc[idx]


def cmdFigures(cfg, args):
    from . import figures as F
    from .experiment import loadBank, makeDataset
    from .features import f0Names
    from .generator import Condition, Generator

    res = resultsDir(cfg)
    figdir = resolve(cfg, "figures")
    bank = loadBank(res / f"bank-{args.source}", cfg)
    mix = tuple(cfg["grid"]["compositions"][-1])
    ds = makeDataset(bank, Condition(c=0.5, classes=mix), cfg, 12345)
    syn = ds["records"][0]
    if args.source == "utd":
        reals = utdAligned(cfg)
    else:
        from .data import loadAww

        reals = loadAww(resolve(cfg, "zhang"))
    real = next(r for r in reals if r.pid == syn.meta["physiologyDonor"])
    F.f1RealVsSynthetic(real, syn, figdir)
    # F2: each class on the same clean segment
    gen = Generator(bank["library"], cfg)
    rng = np.random.default_rng(7)
    m = bank["motion"][0]
    p = next(x for x in bank["phys"] if x["validDonor"] and x["pid"] != m["pid"])
    clean, unused = gen.clean(m, p, rng, 0.0)
    fs = m["fs"]
    seg = slice(int(60 * fs), int(120 * fs))
    ex = {}
    for cls in ("C1", "C2", "C3", "C4", "C5"):
        try:
            prm, unused = gen.drawTemplate(cls, rng) if cls != "C4" else ({}, "")
        except ValueError:
            continue
        n = clean.size
        if cls == "C5":
            dirty = clean.copy()
            a = int(80 * fs)
            dirty[a : a + max(1, int(prm["D"] * fs))] = 0.0
        else:
            ev = {"onset": int(80 * fs), "end": int(90 * fs), "intensity": 1.0}
            if cls == "C4" and not np.isfinite(gen.c4Gain):
                continue
            dirty = clean + gen.artifactSignal(cls, prm, 80.0, 1.0, n, fs, clean, m["acc"], ev, "V1")[0]
        ex[cls] = (clean[seg], dirty[seg])
    F.f2Classes(ex, fs, figdir)
    vj = res / f"validation-{args.source}.json"
    if vj.exists():
        v = json.loads(vj.read_text())
        if "G3" in v:
            F.f3Realism(v["G3"], f0Names, figdir)
    sj = res / f"grid-summary-{args.source}.csv"
    gm = res / f"grid-metrics-{args.source}.csv"
    if sj.exists():
        s = pd.read_csv(sj)
        F.f4Map(s, figdir)  # one file per composition; the mixture is the main map
    if gm.exists():
        df = readGrid(gm, cfg)
        F.f5ByClass(df, figdir)
        F.f6Cost(df, figdir)
    dj = res / "detectors-by-class.csv"
    if dj.exists():
        d = pd.read_csv(dj)
        F.f7Detectors(d.groupby(["detector", "class"]).sensitivity.mean().unstack(), figdir)
    log(f"figures in {figdir}")


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(prog="edasynth", description=__doc__)
    ap.add_argument("--config", default=None, help="YAML config (default configs/default.yaml)")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("check-data")
    s.add_argument("--verbose", action="store_true")
    sub.add_parser("inspect")
    s = sub.add_parser("align")
    s.add_argument("--plots", action="store_true")
    s = sub.add_parser("realdata")
    s.add_argument("--stored", action="store_true")
    s.add_argument("--raw", action="store_true")
    s.add_argument("--out", default=None)
    for name in ("estimate", "detectors", "validate", "pilot", "grid", "analyze", "power", "figures"):
        s = sub.add_parser(name)
        s.add_argument("--source", choices=["utd", "aww"], default="utd")
        s.add_argument("--jobs", type=int, default=1)
        if name == "detectors":
            s.add_argument("--reps", type=int, default=3)
        if name == "validate":
            s.add_argument("--gates", nargs="+", default=["G2", "G3", "G4", "G5", "G6", "G7"])
            s.add_argument("--pytest", action="store_true")
        if name == "grid":
            s.add_argument("--reps", type=int, default=None)
            s.add_argument("--variants", action="store_true")
            s.add_argument("--label-noise", dest="labelNoise", action="store_true")
            s.add_argument("--models", default="logreg")
    args = ap.parse_args(argv)
    cfg = loadConfig(args.config)
    cfg["repoRoot"] = str(Path(__file__).resolve().parents[2])
    fn = globals()["cmd" + "".join(part.capitalize() for part in args.cmd.split("-"))]
    return fn(cfg, args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
