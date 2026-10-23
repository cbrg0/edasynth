"""Figures F1-F7, built only from files in ``results/`` (no manual editing)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CB = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def _save(fig, out: Path, name: str):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def f1_real_vs_synthetic(real, synth, out: Path, minutes: float = 10.0):
    """Real and synthetic recording side by side, same scale (EDA, |a|, labels).

    ``real`` should be the synthetic participant's physiology donor, so that the
    EDA levels are comparable; the synthetic movement comes from another person
    (its motion donor). Each column starts 1 min before its first artifact window.
    """
    fig, ax = plt.subplots(3, 2, figsize=(11, 6), sharex="col", layout="constrained")
    cols = [
        (f"Real ({real.pid})", real.eda, real.acc, real.majority(), real.fs),
        (
            f"Synthetic ({synth.pid}): physiology {synth.meta['physiology_donor']}, "
            f"movement {synth.meta['motion_donor']}",
            synth.y,
            synth.acc,
            synth.labels,
            synth.fs,
        ),
    ]
    for col, (title, eda, acc, lab, fs) in enumerate(cols):
        first = int(np.argmax(lab)) if lab.any() else 0
        a0 = max(0, first * 40 - int(60 * fs))
        n = int(minutes * 60 * fs)
        a1 = min(eda.size, a0 + n)
        t = (np.arange(a0, a1)) / fs / 60
        ax[0, col].plot(t, eda[a0:a1], color=CB[0], lw=0.8)
        ax[0, col].set_title(title, fontsize=9)
        ax[0, col].set_ylabel("EDA, uS")
        ax[1, col].plot(t, np.sqrt((acc[a0:a1] ** 2).sum(1)), color=CB[1], lw=0.6)
        ax[1, col].set_ylabel("|a|, g")
        w = np.repeat(lab, 40)[a0:a1]
        ax[2, col].fill_between(t, 0, w, step="pre", color=CB[3], alpha=0.6)
        ax[2, col].set_ylim(-0.05, 1.05)
        ax[2, col].set_ylabel("artifact")
        ax[2, col].set_xlabel("time, min")
    lo = min(ax[0, 0].get_ylim()[0], ax[0, 1].get_ylim()[0])
    hi = max(ax[0, 0].get_ylim()[1], ax[0, 1].get_ylim()[1])
    for c in (0, 1):
        ax[0, c].set_ylim(lo, hi)
    _save(fig, out, "F1_real_vs_synthetic")


def f2_classes(examples: dict, fs: float, out: Path):
    """Each artifact class added to the same clean segment."""
    fig, ax = plt.subplots(1, len(examples), figsize=(3 * len(examples), 2.6), sharey=True)
    for a, (cls, (clean, dirty)) in zip(np.atleast_1d(ax), examples.items(), strict=True):
        t = np.arange(clean.size) / fs
        a.plot(t, clean, color="0.6", lw=1, label="clean")
        a.plot(t, dirty, color=CB[1], lw=1, label="with artifact")
        a.set_title(cls)
        a.set_xlabel("time, s")
    np.atleast_1d(ax)[0].set_ylabel("EDA, uS")
    np.atleast_1d(ax)[0].legend(frameon=False, fontsize=8)
    _save(fig, out, "F2_artifact_classes")


def f3_realism(g3: dict, names: list[str], out: Path):
    D = np.asarray(g3["ks_D"])
    order = np.argsort(D)[::-1]
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.bar(range(D.size), D[order], color=CB[0])
    ax.set_xticks(range(D.size))
    ax.set_xticklabels([names[k] for k in order], rotation=90, fontsize=7)
    ax.set_ylabel("Kolmogorov-Smirnov D")
    ax.set_title(f"Real vs synthetic clean windows (classifier AUC = {g3['auc']:.2f})")
    _save(fig, out, "F3_realism")


def f4_map(summary: pd.DataFrame, out: Path, reference: dict | None = None):
    """Heat maps of mean dA4 and dB over c x kappa.

    One file per artifact composition (F4_map_<composition>); within a file one
    panel row per (delta_max, pi_r). Compositions are never averaged together.
    """
    s = summary.reset_index()
    vmax = float(np.nanmax(np.abs(s[["mean_dA4", "mean_dB"]].to_numpy()))) or 0.01
    for comp, sc in s.groupby("classes"):
        combos = sc[["delta_max", "pi_r"]].drop_duplicates().to_numpy()
        fig, ax = plt.subplots(
            len(combos), 2, figsize=(10, 3.2 * len(combos) + 0.6), squeeze=False, layout="constrained"
        )
        for i, (dm, pr) in enumerate(combos):
            sub = sc[(sc.delta_max == dm) & (sc.pi_r == pr)]
            for j, d in enumerate(("dA4", "dB")):
                piv = sub.pivot_table(index="kappa", columns="c", values=f"mean_{d}")
                use = sub.pivot_table(index="kappa", columns="c", values=f"useful_{d}")
                a = ax[i, j]
                im = a.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto", origin="lower")
                for (r, c), val in np.ndenumerate(piv.to_numpy()):
                    mark = "*" if bool(use.to_numpy()[r, c]) else ""
                    a.text(c, r, f"{val:+.3f}{mark}", ha="center", va="center", fontsize=8)
                a.set_xticks(range(piv.shape[1]), [f"{v:g}" for v in piv.columns])
                a.set_yticks(range(piv.shape[0]), [f"{v:g}" for v in piv.index])
                a.set_xlabel("c, share of artifacts caused by movement")
                a.set_ylabel("kappa")
                a.set_title(f"{'F0+FA4' if d == 'dA4' else 'F0+FB'}; dmax = {dm:g} s, pi_r = {pr:g}", fontsize=9)
                fig.colorbar(im, ax=a, label="mean delta AUPRC", shrink=0.9)
                if reference is not None and np.isfinite(reference.get("c", np.nan)):
                    a.axvline(np.interp(reference["c"], piv.columns, range(piv.shape[1])), color="k", ls="--", lw=1)
        note = "; dashed line: UTD label reference (biased)" if reference is not None else ""
        fig.suptitle(f"Artifacts {comp}. * = mean >= 0.01 and lower 95 % bound > 0{note}", fontsize=10)
        _save(fig, out, f"F4_map_{comp.replace('+', '-')}")


def f5_by_class(df: pd.DataFrame, out: Path):
    single = df[~df.classes.str.contains(r"\+")]
    fig, ax = plt.subplots(1, single.classes.nunique(), figsize=(3.2 * single.classes.nunique(), 2.8), sharey=True)
    for a, (cls, sub) in zip(np.atleast_1d(ax), single.groupby("classes"), strict=False):
        for k, d in enumerate(("dA4", "dB")):
            g = sub.groupby("c")[d]
            m, s, n = g.mean(), g.std(), g.count()
            a.errorbar(m.index, m, yerr=1.96 * s / np.sqrt(n), color=CB[k], marker="o", label=d, capsize=2)
        a.axhline(0, color="k", lw=0.8)
        a.axhline(0.01, color="0.5", ls=":", lw=0.8)
        a.set_title(cls)
        a.set_xlabel("c")
    np.atleast_1d(ax)[0].set_ylabel("delta AUPRC")
    np.atleast_1d(ax)[0].legend(frameon=False)
    _save(fig, out, "F5_by_class")


def f6_cost(df: pd.DataFrame, out: Path):
    sub = df[df.pi_r > 0]
    cols = ["lost_motion_F0", "lost_motion_F0+FA4", "lost_motion_F0+FB"]
    m = sub[cols].mean()
    s = sub[cols].std() / np.sqrt(len(sub)) * 1.96
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.bar(range(3), m, yerr=s, color=CB[:3], capsize=3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["F0", "F0+FA4", "F0+FB"])
    ax.set_ylabel("share of real motion responses\nflagged as artifact")
    ax.set_title("Cost of the accelerometer (pi_r = 0.5)")
    _save(fig, out, "F6_cost")


def f7_detectors(table: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(6, 3))
    im = ax.imshow(table.to_numpy(float), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    for (r, c), v in np.ndenumerate(table.to_numpy(float)):
        ax.text(c, r, f"{v:.2f}", ha="center", va="center", color="w" if v < 0.6 else "k", fontsize=8)
    ax.set_xticks(range(table.shape[1]))
    ax.set_xticklabels(table.columns)
    ax.set_yticks(range(table.shape[0]))
    ax.set_yticklabels(table.index)
    fig.colorbar(im, ax=ax, label="share of artifact windows flagged")
    _save(fig, out, "F7_detectors_by_class")
