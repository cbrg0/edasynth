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
from .estimation import decompose_participant, response_rates, scr_events
from .experiment import evaluate, make_dataset
from .features import f0
from .generator import Condition, Generator, derangement


# ------------------------------------------------------------------ G2
def g2_parameter_recovery(bank: dict, cfg: dict, n: int, seed: int) -> dict:
    """Generate artifact-free participants with known physiology; re-estimate it."""
    v = cfg["validation"]
    rng = np.random.default_rng(seed)
    gen = Generator(bank["library"], cfg)
    cond = Condition(prevalence=0.0, pi_r=0.0)
    m_all, p_all = bank["motion"], bank["phys"]
    perm = derangement(len(m_all), rng)
    tau_err, rate_rel, true_a, est_a = [], [], [], []
    for k in range(min(n, len(m_all))):
        m, p = m_all[k], p_all[int(perm[k])]
        rec = gen.generate(f"g2_{k}", m, p, cond, rng)
        nw = rec.labels.size
        r = Recording(pid=rec.pid, fs=rec.fs, eda=rec.y, acc=rec.acc, stage=rec.stage, labels3=np.zeros((nw, 3), int))
        clean = np.ones(nw, bool)
        dec = decompose_participant(r, clean, cfg)
        ev = scr_events(
            dec["r"], r.fs, dec["tau0"], cfg["estimation"]["cvxeda"]["tau1"], cfg["estimation"]["min_response_amp"]
        )
        rates = response_rates(r, ev, clean)
        tau_err.append(abs(dec["tau0"] - p["tau0"]))
        for g in ("rest", "load"):
            t = sum(
                1
                for x in rec.responses
                if x["origin"] == "spontaneous" and x["amp"] >= 0.05 and _stage_type_at(rec, x["t"]) == g
            )
            T = rates[g]["T_s"]
            if T > 0 and t > 0:
                rate_rel.append(rates[g]["rate_per_s"] / (t / T) - 1)
        true_a += [x["area"] for x in rec.responses if x["amp"] >= 0.05]
        est_a += [e["area"] for e in ev]
    ks = stats.ks_2samp(true_a, est_a) if true_a and est_a else None
    med_tau = float(np.median(tau_err))
    med_rate = float(np.median(np.abs(rate_rel))) if rate_rel else float("nan")
    ok = (
        med_tau <= v["g2_tau0_max_err_s"]
        and med_rate <= v["g2_rate_max_rel"]
        and ks is not None
        and ks.pvalue > v["g2_ks_p"]
    )
    return {
        "passed": bool(ok),
        "median_tau0_error_s": med_tau,
        "median_abs_rate_rel_error": med_rate,
        "ks_area_p": float(ks.pvalue) if ks else float("nan"),
        "ks_area_D": float(ks.statistic) if ks else float("nan"),
        "n_participants": len(tau_err),
    }


def _stage_type_at(rec, t):
    from .data import stage_type

    i = min(int(t * rec.fs), rec.stage.size - 1)
    return stage_type(str(rec.stage[i]))


# ------------------------------------------------------------------ G3, G4
def classifier_two_sample(Xr, gr, Xs, gs, folds: int = 5) -> dict:
    """AUC of a classifier separating real from synthetic windows (grouped CV)."""
    X = np.vstack([Xr, Xs])
    y = np.r_[np.zeros(len(Xr)), np.ones(len(Xs))]
    g = np.r_[gr, gs + 10_000]
    pred = np.zeros(y.size)
    for tr, te in GroupKFold(n_splits=folds).split(X, y, g):
        m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)).fit(X[tr], y[tr])
        pred[te] = m.predict_proba(X[te])[:, 1]
    D = [float(stats.ks_2samp(Xr[:, k], Xs[:, k]).statistic) for k in range(X.shape[1])]
    return {"auc": float(roc_auc_score(y, pred)), "ks_D": D}


def g3_g4_realism(real: list[Recording], bank: dict, cfg: dict, seed: int) -> dict:
    v = cfg["validation"]
    deriv = cfg["features"]["derivative"]
    mix = tuple(cfg["grid"]["compositions"][-1])
    ds = make_dataset(bank, Condition(c=0.5, classes=mix, prevalence=cfg["grid"]["prevalence"]), cfg, seed)
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
    g3 = classifier_two_sample(Xr[unan], gr[unan], Xs[ys == 0], gs[ys == 0])
    g4 = classifier_two_sample(Xr[yr == 1], gr[yr == 1], Xs[ys == 1], gs[ys == 1])
    g3["passed"] = g3["auc"] <= v["c2st_auc_max"]
    g4["passed"] = g4["auc"] <= v["c2st_auc_max"]
    g4["note"] = "UTD artifact labels partly depend on movement; interpret with that caveat"
    return {"G3": g3, "G4": g4}


# ------------------------------------------------------------------ G5, G6
def _reps(bank, cond, cfg, reps, seed0, acc_shift=False):
    rows = []
    for r in range(reps):
        ds = make_dataset(bank, cond, cfg, seed0 + r, acc_shift=acc_shift)
        row = evaluate(ds, cfg, ("logreg",), seed0 + r)[0]
        row.pop("_preds")
        rows.append(row)
    return rows


def _t_ci(x):
    x = np.asarray(x, float)
    m = x.mean()
    h = stats.t.ppf(0.975, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size)
    return float(m), float(m - h), float(m + h)


def g5_negative_controls(bank, cfg, seed) -> dict:
    v = cfg["validation"]
    mix = tuple(cfg["grid"]["compositions"][-1])
    out = {}
    ok = True
    for name, cond, shift in (
        ("c0", Condition(c=0.0, classes=mix), False),
        ("acc_shift", Condition(c=1.0, kappa=1.0, delta_max=1.0, classes=mix), True),
    ):
        rows = _reps(bank, cond, cfg, v["control_reps"], seed, acc_shift=shift)
        res = {}
        for d in ("dA", "dA4", "dB"):
            m, lo, hi = _t_ci([r[d] for r in rows])
            good = lo <= 0 <= hi and abs(m) < v["g5_max_abs_delta"]
            res[d] = {"mean": m, "ci": [lo, hi], "passed": bool(good)}
            ok &= good
        out[name] = res
    out["passed"] = bool(ok)
    return out


def g6_positive_controls(bank, cfg, seed) -> dict:
    v = cfg["validation"]
    out = {}
    ok = True
    for name, cond in (
        ("c4_only", Condition(c=1.0, classes=("C4",))),
        (
            "c_1_kappa_1_dmax_1_mixture",
            Condition(c=1.0, kappa=1.0, delta_max=1.0, classes=tuple(cfg["grid"]["compositions"][-1])),
        ),
    ):
        rows = _reps(bank, cond, cfg, v["control_reps"], seed)
        res = {"n_positive_windows": [r["n_positive"] for r in rows]}
        any_ok = False
        if min(res["n_positive_windows"]) == 0:
            res["note"] = "no positive windows were generated; the control cannot be evaluated"
        for d in ("dA", "dB"):
            m, lo, hi = _t_ci([r[d] for r in rows])
            res[d] = {"mean": m, "ci": [lo, hi]}
            any_ok |= m > 0 and lo > 0
        res["passed"] = bool(any_ok)
        ok &= any_ok
        out[name] = res
    out["passed"] = bool(ok)
    return out


def g7_real_consistency(raw_result: dict, cfg: dict) -> dict:
    v = cfg["validation"]
    ref = v["g7_reference"]
    same_sign = np.sign(raw_result["delta"]) == np.sign(ref["delta"])
    close = (
        abs(raw_result["auprc_0"] - ref["auprc_0"]) <= v["g7_max_abs_diff"]
        and abs(raw_result["auprc_1"] - ref["auprc_1"]) <= v["g7_max_abs_diff"]
    )
    return {
        "passed": bool(same_sign and close),
        "ours": {k: raw_result[k] for k in ("auprc_0", "auprc_1", "delta")},
        "reference": ref,
    }
