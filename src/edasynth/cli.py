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
from .config import load_config, resolve

RESULTS_DEFAULT = "results"


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _write_json(path: Path, obj):
    from .bank import to_json

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(to_json(obj))


def _results(cfg) -> Path:
    p = resolve(cfg, "results")
    p.mkdir(parents=True, exist_ok=True)
    return p


# ------------------------------------------------------------------ data
def cmd_check_data(cfg, args):
    """Verify SHA-256 checksums of the downloaded files."""
    bad = 0
    lists = (
        ("noneeg", resolve(cfg, "noneeg") / "SHA256SUMS.txt"),  # official list shipped with the data
        ("zhang", Path(cfg["_repo"]) / "data/checksums/zhang-eda-artifact-detection-7c4bd46.sha256"),
    )
    for key, listfile in lists:
        root = resolve(cfg, key)
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
                print(f"MISSING  {key}/{name}")
                bad += 1
                continue
            got = hashlib.sha256(f.read_bytes()).hexdigest()
            ok = got == digest
            bad += not ok
            if not ok or args.verbose:
                print(f"{'OK      ' if ok else 'MISMATCH'} {key}/{name}")
    print("all files verified" if bad == 0 else f"{bad} problem(s)")
    return 1 if bad else 0


def _utd(cfg, with_labels=True):
    from .data import load_utd_labels, load_utd_participant

    labels = load_utd_labels(resolve(cfg, "zhang")) if with_labels else None
    excl = set(cfg["data"]["exclude_subjects"])
    recs = []
    for s in range(1, 21):
        if s in excl:
            continue
        recs.append(load_utd_participant(resolve(cfg, "noneeg"), s, labels))
    return recs


def _utd_aligned(cfg):
    """UTD recordings accepted by `align`, with the alignment offset applied."""
    path = _results(cfg) / "alignment.csv"
    if not path.exists():
        raise SystemExit("results/alignment.csv not found: run `python -m edasynth align` first")
    align = pd.read_csv(path)
    ok = align[align.accepted]
    off = dict(zip(ok.pid, ok.offset_samples, strict=True))
    return [r.cropped(int(off[r.pid])) for r in _utd(cfg) if r.pid in off]


def cmd_inspect(cfg, args):
    """Units, ranges, annotation codes and stage parsing for every UTD record."""
    from .data import stage_type

    rows = []
    for r in _utd(cfg):
        mag = np.sqrt((r.acc**2).sum(1))
        rest = np.array([stage_type(str(s)) == "rest" for s in r.stage]) if r.stage is not None else np.ones(r.n, bool)
        rows.append(
            {
                "pid": r.pid,
                "n_samples": r.n,
                "n_windows": r.n_windows,
                "fs": r.fs,
                "units": r.units,
                "eda_min": float(r.eda.min()),
                "eda_median": float(np.median(r.eda)),
                "eda_max": float(r.eda.max()),
                "acc_mag_median_rest": float(np.median(mag[rest])),
                "n_annotations": len(r.meta.get("annotations", [])),
                "annotation_codes": sorted({a[1] for a in r.meta.get("annotations", [])}),
                "annotation_aux": sorted({a[2] for a in r.meta.get("annotations", [])}),
                "stage_error": r.meta.get("stage_error", ""),
            }
        )
    out = _results(cfg) / "inspect_utd.json"
    _write_json(out, rows)
    df = pd.DataFrame(rows)
    print(
        df[
            [
                "pid",
                "n_samples",
                "eda_min",
                "eda_median",
                "eda_max",
                "acc_mag_median_rest",
                "n_annotations",
                "stage_error",
            ]
        ].to_string()
    )
    log(f"written {out}")


def cmd_align(cfg, args):
    """Check that the expert labels line up with the raw signals (docs/METHOD.md section 2.4)."""
    from .alignment import check_participant
    from .data import load_utd_labels, load_utd_stored_features

    labels = load_utd_labels(resolve(cfg, "zhang"))
    eda_stored, _ = load_utd_stored_features(resolve(cfg, "zhang"))
    rows = []
    for r in _utd(cfg, with_labels=False):
        g = int(r.pid[1:]) - 1
        m = labels.group.to_numpy() == g
        res = check_participant(r.eda, eda_stored[m][:, :4], cfg["data"]["alignment_min_r"])
        res.update(pid=r.pid, group=g, artifact_share=float(labels.majority[m].mean()))
        rows.append(res)
        log(
            f"{r.pid}: windows raw/labels {res['n_windows_raw']}/{res['n_windows_labels']}, "
            f"r_min={res['r_min']:.4f}, shift={res['shift_windows']}w+{res['shift_samples']}s, "
            f"{'ACCEPTED' if res['accepted'] else 'REJECTED'}"
        )
    df = pd.DataFrame(rows)
    out = _results(cfg) / "alignment.csv"
    df.to_csv(out, index=False)
    if args.plots:
        import matplotlib.pyplot as plt

        from .figures import CB, _save

        offsets = dict(zip(df.pid, df.offset_samples, strict=True))
        for r0 in _utd(cfg):
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
            _save(fig, resolve(cfg, "figures") / "alignment", f"step2_{r.pid}")
    log(f"written {out}; accepted {int(df.accepted.sum())}/{len(df)}")


# ------------------------------------------------------------------ real data
def cmd_realdata(cfg, args):
    """Baseline on real UTD: stored features (reproduction) and/or raw recomputation (G7)."""
    from .data import load_utd_labels, load_utd_stored_features
    from .features import f0, fa4
    from .realdata import compare, shuffled_control, single_feature_auprc

    res = {}
    labels = load_utd_labels(resolve(cfg, "zhang"))
    y = labels.majority.to_numpy()
    g = labels.group.to_numpy()
    rcfg = _realdata_cfg(cfg)
    if args.stored:
        eda, acc = load_utd_stored_features(resolve(cfg, "zhang"))
        log("stored features: logistic regression, LOSO")
        res["stored_logreg"] = compare(eda, np.hstack([eda, acc[:, :4]]), y, g, rcfg, "logreg")
        log("stored features: random forest, LOSO")
        res["stored_rf"] = compare(eda, np.hstack([eda, acc[:, :4]]), y, g, rcfg, "rf")
        log("stored features: shuffled-acceleration control")
        res["stored_shuffled"] = shuffled_control(eda, acc[:, :4], y, g, rcfg, cfg["realdata"]["n_shuffles"], 0)
        res["stored_single_feature_auprc"] = single_feature_auprc(eda, y)
    if args.raw:
        from scipy import stats

        acc_rows, eda_rows, yy, gg = [], [], [], []
        for r in _utd_aligned(cfg):
            nw = r.labels3.shape[0]
            # z-scored within participant, as in the published feature matrices
            eda_rows.append(np.nan_to_num(stats.zscore(f0(r.eda, cfg["features"]["derivative"])[:nw], axis=0)))
            acc_rows.append(np.nan_to_num(stats.zscore(fa4(r.acc)[:nw], axis=0)))
            yy.append(r.majority())
            gg.append(np.full(nw, int(r.pid[1:])))
        E, A = np.vstack(eda_rows), np.vstack(acc_rows)
        Y, G = np.concatenate(yy), np.concatenate(gg)
        res["raw_logreg"] = compare(E, np.hstack([E, A]), Y, G, rcfg, "logreg")
        from .validation import g7_real_consistency

        res["G7"] = g7_real_consistency(res["raw_logreg"], cfg)
    for v in res.values():
        if isinstance(v, dict):
            v.pop("_p0", None)
            v.pop("_p1", None)
    out = _results(cfg) / ("realdata.json" if not args.out else args.out)
    _write_json(out, res)
    for k, v in res.items():
        if isinstance(v, dict) and "delta" in v:
            log(
                f"{k}: AUPRC {v['auprc_0']:.4f} -> {v['auprc_1']:.4f}, delta {v['delta']:+.4f}, "
                f"95% CI [{v['ci95'][0]:+.4f}, {v['ci95'][1]:+.4f}], Wilcoxon p={v['wilcoxon_p']:.3f}"
            )
    log(f"written {out}")


def _realdata_cfg(cfg):
    from .config import merge

    return merge(cfg, {"models": cfg["realdata"]["models"]})


# ------------------------------------------------------------------ estimation
def cmd_estimate(cfg, args):
    """Estimate generator parameters (docs/METHOD.md section 3) and write the bank."""
    from .bank import build_bank

    if args.source == "utd":
        recs = _utd_aligned(cfg)
        bad = [r.pid for r in recs if r.stage is None]
        if bad:
            raise SystemExit(f"stages could not be parsed for {bad}; run `inspect` first")
    else:
        from .data import load_aww

        recs = load_aww(resolve(cfg, "zhang"))
        log("AWW is ONE person: donors are different recording days, not different people (development only)")
    out = _results(cfg) / f"bank_{args.source}"
    meta = build_bank(recs, cfg, args.source, out, log, n_jobs=args.jobs)
    lib = meta["library"]
    log(f"episodes by class: {lib['episode_counts']}; C4 gain G = {lib['G']:.4f}")
    for e in meta.get("estimates", []):
        for w in e.get("warnings", []):
            log(f"WARNING {e['pid']}: {w}")
    log(f"written {out}.json / .npz")


def _bank(cfg, args):
    from .experiment import load_bank

    return load_bank(_results(cfg) / f"bank_{args.source}")


# ------------------------------------------------------------------ detectors (goal 2)
def cmd_detectors(cfg, args):
    from .detectors import run_all
    from .experiment import make_dataset
    from .generator import Condition
    from .models import auprc

    bank = _bank(cfg, args)
    rows = []
    for ci, cls in enumerate(cfg["detectors"]["classes"]):
        for rep in range(args.reps):
            cond = Condition(c=0.5, classes=(cls,), prevalence=cfg["grid"]["prevalence"])
            ds = make_dataset(bank, cond, cfg, cfg["seeds"]["detectors"] + rep * 31 + ci * 1000)
            flags = {}
            for r in ds["records"]:
                for name, (kind, val) in run_all(r.y, r.fs, cfg).items():
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
                        "n_positive": 0,
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
                            "n_positive": int(y.sum()),
                        }
                    )
                else:
                    from .models import sensitivity_at_specificity

                    rows.append(
                        {
                            "class": cls,
                            "rep": rep,
                            "detector": name,
                            "sensitivity": sensitivity_at_specificity(y, s, 0.95),
                            "specificity": 0.95,
                            "auprc": auprc(y, s),
                            "n_positive": int(y.sum()),
                        }
                    )
        log(f"detectors done for {cls}")
    df = pd.DataFrame(rows)
    out = _results(cfg) / "detectors_by_class.csv"
    df.to_csv(out, index=False)
    print(df.groupby(["detector", "class"]).sensitivity.mean().unstack().round(3).to_string())
    log(f"written {out}")


# ------------------------------------------------------------------ validation
def cmd_validate(cfg, args):
    from . import validation as V

    bank = _bank(cfg, args)
    rep = {"source": args.source, "version": __version__}
    if args.pytest:
        r = subprocess.run([sys.executable, "-m", "pytest", "-q"], capture_output=True, text=True)
        rep["G1"] = {"passed": r.returncode == 0, "tail": r.stdout.strip().splitlines()[-1:]}
    seed = cfg["seeds"]["validation"]
    out = _results(cfg) / f"validation_{args.source}.json"

    def save():  # written after every gate, so a crash later does not lose finished gates
        _write_json(out, rep)

    if "G2" in args.gates:
        log("G2 parameter recovery")
        rep["G2"] = V.g2_parameter_recovery(bank, cfg, cfg["validation"]["g2_participants"], seed)
        save()
    if "G3" in args.gates or "G4" in args.gates:
        if args.source == "utd":
            log("G3/G4 realism")
            rep.update(V.g3_g4_realism(_utd_aligned(cfg), bank, cfg, seed))
        else:
            from .data import load_aww

            log("G3/G4 realism against AWW (development only)")
            rep.update(V.g3_g4_realism(load_aww(resolve(cfg, "zhang")), bank, cfg, seed))
        save()
    if "G5" in args.gates:
        log("G5 negative controls")
        rep["G5"] = V.g5_negative_controls(bank, cfg, seed + 1)
        save()
    if "G6" in args.gates:
        log("G6 positive controls")
        rep["G6"] = V.g6_positive_controls(bank, cfg, seed + 2)
        save()
    rj = _results(cfg) / "realdata.json"
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
def _run_jobs(cfg, bank, jobs, n_jobs, out_csv, models, label_noise=None):
    from joblib import Parallel, delayed

    from .experiment import run_cell

    done = set()
    if out_csv.exists():
        old = pd.read_csv(out_csv)
        done = set(zip(old["key"], old["rep"], strict=True))
    todo = [(c, r) for c, r in jobs if (c.key(), r) not in done]
    log(f"{len(todo)} replicate-cells to run ({len(jobs) - len(todo)} already done)")
    batch = max(1, n_jobs * 2)
    for i in range(0, len(todo), batch):
        chunk = todo[i : i + batch]
        res = Parallel(n_jobs=n_jobs)(delayed(run_cell)(bank, c, cfg, r, models, label_noise) for c, r in chunk)
        rows = []
        for (c, _r), rr in zip(chunk, res, strict=True):
            for row in rr:
                row["key"] = c.key()
                rows.append(row)
        df = pd.DataFrame(rows)
        # rows of different cells can have different columns (e.g. det_C2_* only for
        # mixtures), so the file is rewritten with aligned columns instead of appended
        if out_csv.exists():
            df = pd.concat([pd.read_csv(out_csv), df], ignore_index=True, sort=False)
        tmp = out_csv.with_suffix(".tmp")
        df.to_csv(tmp, index=False)
        tmp.replace(out_csv)
        log(f"{min(i + batch, len(todo))}/{len(todo)}")


def cmd_pilot(cfg, args):
    from .experiment import required_reps
    from .generator import Condition

    bank = _bank(cfg, args)
    mix = tuple(cfg["grid"]["compositions"][-1])
    conds = [Condition(c=c, kappa=k, classes=mix) for c in (0.0, 1.0) for k in (0.0, 1.0)]
    jobs = [(c, r) for c in conds for r in range(cfg["pilot"]["reps"])]
    out = _results(cfg) / f"pilot_{args.source}.csv"
    _run_jobs(cfg, bank, jobs, args.jobs, out, ("logreg",))
    df = pd.read_csv(out)
    sd = float(df.groupby("key")[["dA4", "dB"]].std().max().max())
    R = required_reps(sd, cfg["pilot"]["half_width"], cfg["pilot"]["max_reps"])
    _write_json(_results(cfg) / f"pilot_{args.source}_R.json", {"max_sd_delta": sd, "R": R})
    log(f"largest SD of delta across pilot cells = {sd:.4f} -> R = {R}")


def cmd_grid(cfg, args):
    from .experiment import grid_conditions, variant_conditions

    bank = _bank(cfg, args)
    R = args.reps
    if R is None:
        R = json.loads((_results(cfg) / f"pilot_{args.source}_R.json").read_text())["R"]
    conds = variant_conditions(cfg) if args.variants else grid_conditions(cfg)
    models = tuple(args.models.split(","))
    jobs = [(c, r) for c in conds for r in range(R)]
    name = "variants" if args.variants else "grid"
    if args.label_noise:
        name += "_labelnoise"
    out = _results(cfg) / f"{name}_metrics_{args.source}.csv"
    _run_jobs(
        cfg, bank, jobs, args.jobs, out, models, cfg["analysis"]["label_noise_kappa"] if args.label_noise else None
    )
    log(f"written {out}")


def cmd_analyze(cfg, args):
    from .experiment import compatibility, summarise

    res = _results(cfg)
    df = pd.read_csv(res / f"grid_metrics_{args.source}.csv")
    df = df[df.model == "logreg"]
    by = ["c", "kappa", "delta_max", "pi_r", "classes"]
    parts = []
    for d in ("dA", "dA4", "dB", "dB_minus_dA4"):
        s = summarise(df, by, d, cfg["analysis"]["decision"])
        s.columns = [f"{c}_{d}" for c in s.columns]
        parts.append(s)
    summ = pd.concat(parts, axis=1)
    summ.to_csv(res / f"grid_summary_{args.source}.csv")
    comp = compatibility(
        df[df.classes == "+".join(cfg["grid"]["compositions"][-1])],
        cfg["analysis"]["observed_delta"],
        ["c", "kappa", "delta_max", "pi_r"],
    )
    comp.to_csv(res / f"compatibility_{args.source}.csv")
    ref = cfg["analysis"]["sensitivity_reference"]
    mix = "+".join(cfg["grid"]["compositions"][-1])
    sens = []
    base = df[df.classes == mix]
    for p, levels in (
        ("c", cfg["grid"]["c"]),
        ("kappa", cfg["grid"]["kappa"]),
        ("delta_max", cfg["grid"]["delta_max"]),
        ("pi_r", cfg["grid"]["pi_r"]),
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
                "range_dA4": float(means.dA4.max() - means.dA4.min()),
                "range_dB": float(means.dB.max() - means.dB.min()),
            }
        )
    pd.DataFrame(sens).sort_values("range_dB", ascending=False).to_csv(
        res / f"sensitivity_{args.source}.csv", index=False
    )
    log(f"compatible cells: {int(comp.compatible.sum())}/{len(comp)}")
    print(pd.DataFrame(sens).to_string())


def cmd_power(cfg, args):
    from .experiment import bootstrap_power
    from .generator import Condition

    bank = _bank(cfg, args)
    res = _results(cfg)
    comp = pd.read_csv(res / f"compatibility_{args.source}.csv")
    comp = comp[comp.compatible].head(cfg["power"]["n_conditions"])
    mix = tuple(cfg["grid"]["compositions"][-1])
    rows = []
    for _, c in comp.iterrows():
        cond = Condition(c=c.c, kappa=c.kappa, delta_max=c.delta_max, pi_r=c.pi_r, classes=mix)
        for n in cfg["power"]["participants"]:
            pw = bootstrap_power(bank, cond, cfg, n, cfg["power"]["reps"])
            rows.append({"condition": cond.key(), "n_participants": n, **{f"power_{k}": v for k, v in pw.items()}})
            log(f"{cond.key()} n={n}: {pw}")
    pd.DataFrame(rows).to_csv(res / f"protocol_planning_{args.source}.csv", index=False)


def cmd_figures(cfg, args):
    from . import figures as F
    from .experiment import load_bank, make_dataset
    from .features import F0_NAMES
    from .generator import Condition, Generator

    res = _results(cfg)
    figdir = resolve(cfg, "figures")
    bank = load_bank(res / f"bank_{args.source}")
    mix = tuple(cfg["grid"]["compositions"][-1])
    ds = make_dataset(bank, Condition(c=0.5, classes=mix), cfg, 12345)
    syn = ds["records"][0]
    if args.source == "utd":
        reals = _utd_aligned(cfg)
    else:
        from .data import load_aww

        reals = load_aww(resolve(cfg, "zhang"))
    real = next(r for r in reals if r.pid == syn.meta["physiology_donor"])
    F.f1_real_vs_synthetic(real, syn, figdir)
    # F2: each class on the same clean segment
    gen = Generator(bank["library"], cfg)
    rng = np.random.default_rng(7)
    m, p = bank["motion"][0], bank["phys"][1 % len(bank["phys"])]
    clean, _ = gen.clean(m, p, rng, 0.0)
    fs = m["fs"]
    seg = slice(int(60 * fs), int(120 * fs))
    ex = {}
    for cls in ("C1", "C2", "C3", "C4", "C5"):
        try:
            prm, _ = gen._template(cls, rng) if cls != "C4" else ({}, "")
        except ValueError:
            continue
        n = clean.size
        if cls == "C5":
            dirty = clean.copy()
            a = int(80 * fs)
            dirty[a : a + max(1, int(prm["D"] * fs))] = 0.0
        else:
            ev = {"onset": int(80 * fs), "end": int(90 * fs), "intensity": 1.0}
            if cls == "C4" and not np.isfinite(bank["library"]["G"]):
                continue
            dirty = clean + gen._artifact_signal(cls, prm, 80.0, 1.0, n, fs, clean, m["acc"], ev, "V1")[0]
        ex[cls] = (clean[seg], dirty[seg])
    F.f2_classes(ex, fs, figdir)
    vj = res / f"validation_{args.source}.json"
    if vj.exists():
        v = json.loads(vj.read_text())
        if "G3" in v:
            F.f3_realism(v["G3"], F0_NAMES, figdir)
    sj = res / f"grid_summary_{args.source}.csv"
    gm = res / f"grid_metrics_{args.source}.csv"
    if sj.exists():
        s = pd.read_csv(sj)
        F.f4_map(s, figdir)  # one file per composition; the mixture is the main map
    if gm.exists():
        df = pd.read_csv(gm)
        df = df[df.model == "logreg"]
        F.f5_by_class(df, figdir)
        F.f6_cost(df, figdir)
    dj = res / "detectors_by_class.csv"
    if dj.exists():
        d = pd.read_csv(dj)
        F.f7_detectors(d.groupby(["detector", "class"]).sensitivity.mean().unstack(), figdir)
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
            s.add_argument("--label-noise", action="store_true")
            s.add_argument("--models", default="logreg")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    cfg["_repo"] = str(Path(__file__).resolve().parents[2])
    fn = globals()["cmd_" + args.cmd.replace("-", "_")]
    return fn(cfg, args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
