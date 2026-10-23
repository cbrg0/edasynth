"""Estimating generator parameters from real labelled recordings (method step 3).

Everything the literature does not provide (noise level, quantisation step,
response rate and size, tonic amplitudes, artifact shapes, motion events) is
measured here on real wrist recordings. The generator later resamples these
measured values (empirical distributions) instead of assuming a parametric form.
"""

from __future__ import annotations

import numpy as np

from .data import WIN, Recording, stage_type
from .decompose import decompose, h_max
from .features import haar_details, jerk, windows


# ---------------------------------------------------------------- clean windows
def clean_windows(rec: Recording) -> np.ndarray:
    """Windows that all three experts called clean AND whose neighbours did too.

    One window on each side corresponds to the 5-s transition zone of
    Kleckner et al. (2018); label errors concentrate at artifact borders.
    """
    unan_clean = rec.labels3.sum(1) == 0
    ok = unan_clean.copy()
    ok[1:] &= unan_clean[:-1]
    ok[:-1] &= unan_clean[1:]
    return ok


def window_mask_to_samples(mask: np.ndarray, n: int) -> np.ndarray:
    m = np.zeros(n, bool)
    m[: mask.size * WIN] = np.repeat(mask, WIN)
    return m


def window_stage(rec: Recording) -> np.ndarray:
    """Stage of each window (stage of its first sample)."""
    return rec.stage[: rec.n_windows * WIN : WIN]


# ---------------------------------------------------------------- noise, quantisation
def noise_sigma(rec: Recording, clean: np.ndarray, with_source: bool = False):
    """sigma = median(|d1|) / 0.6745 over clean rest windows (Donoho & Johnstone, 1994).

    A recording without any clean rest window (e.g. an AWW walking day) falls back
    to all clean windows; the source is reported so the bank records it. Movement
    windows may contain fast physiological or mechanical changes, so the fallback
    can overestimate sigma.
    """
    st = window_stage(rec)
    sel = clean & np.array([stage_type(s) == "rest" for s in st])
    source = "clean_rest"
    if not sel.any():
        sel, source = clean.copy(), "clean_all_fallback"
    w = windows(rec.eda)[sel]
    if w.size == 0:
        sigma, source = float("nan"), "none"
    else:
        d1 = haar_details(w)[0]
        sigma = float(np.median(np.abs(d1)) / 0.6745)
    return (sigma, source) if with_source else sigma


def quantisation_step(x: np.ndarray, decimals: int = 9, min_share: float = 0.01) -> float:
    """Quantisation step of a recorded signal.

    The smallest non-zero absolute difference between neighbouring samples that
    occurs in at least ``min_share`` of all non-zero differences (rare tiny values
    from rounding noise are ignored).
    """
    d = np.round(np.abs(np.diff(x)), decimals)
    d = d[d > 0]
    if d.size == 0:
        return 0.0
    vals, counts = np.unique(d, return_counts=True)
    ok = vals[counts >= min_share * d.size]
    return float(ok.min()) if ok.size else float(vals[np.argmax(counts)])


# ---------------------------------------------------------------- decomposition
def decompose_participant(rec: Recording, clean: np.ndarray, cfg: dict) -> dict:
    """z-score the whole record, decompose with each tau0 on the grid, keep the best.

    tau0 is chosen by the smallest value of the cvxEDA objective (fit plus
    sparsity and smoothness penalties) on the whole record [decision]. The
    residual norm alone was tried first and rejected: it is almost flat in tau0
    and pointed to the wrong end of the grid on synthetic records with known
    tau0 (tests/test_estimation.py::test_tau0_selection_recovers_known_value).
    The residual norm on clean samples is still stored for reference. The
    result is converted back to signal units.
    """
    c = cfg["estimation"]["cvxeda"]
    y = rec.eda.astype(float)
    mu, sd = float(y.mean()), float(y.std())
    z = (y - mu) / sd
    mask = window_mask_to_samples(clean, rec.n)
    best, p_warm = None, None
    grid = np.round(np.arange(c["tau0_min"], c["tau0_max"] + 1e-9, c["tau0_step"]), 3)
    scores = {}
    for tau0 in grid:
        out = decompose(
            z,
            1.0 / rec.fs,
            float(tau0),
            c["tau1"],
            c["alpha"],
            c["gamma"],
            c["knot_s"],
            backend=c["backend"],
            p0=p_warm,
        )
        p_warm = out["p"]
        err = float(out["obj"])
        scores[float(tau0)] = {"objective": err, "residual_clean": float(np.linalg.norm(out["e"][mask]))}
        if best is None or err < best[1]:
            best = (float(tau0), err, out)
    tau0, _, out = best
    return {
        "tau0": tau0,
        "tau0_scores": scores,
        "kkt": out.get("kkt"),
        "iterations": out.get("iterations"),
        "r": sd * out["r"],
        "p": sd * out["p"],
        "t": sd * out["t"] + mu,
        "backend": out["backend"],
        "tau0_at_bound": tau0 in (grid[0], grid[-1]),
    }


def scr_events(r: np.ndarray, fs: float, tau0: float, tau1: float, min_amp: float) -> list[dict]:
    """Responses as peaks of the phasic component r with trough-to-peak height >= ``min_amp``.

    Height A is the peak's prominence (rise above the deepest point between it and
    the neighbouring higher peaks); onset is the latest lowest sample in the 5 s
    before the peak. The area that reproduces this height in the generator is
    a = A / h_max(tau0).

    Counting runs of non-zero driver p was tried first and rejected: on noisy,
    quantised recordings the cvxEDA driver is dense (most samples non-zero at the
    optimum), so runs fragment into many spurious "responses" (AWW development
    data: 14-95 per minute). Peaks of r gave 1.7-2.5 per minute on the AWW
    resting days.
    """
    from scipy.signal import find_peaks

    pk, prop = find_peaks(r, prominence=min_amp)
    hm = h_max(tau0, tau1)
    back = int(round(5 * fs))
    events = []
    for k, A in zip(pk, prop["prominences"], strict=True):
        s0 = max(0, k - back)
        seg = r[s0 : k + 1]
        onset = s0 + seg.size - 1 - int(np.argmin(seg[::-1]))  # last (latest) minimum before the peak
        events.append({"onset": int(onset), "peak": int(k), "area": float(A / hm), "amp": float(A)})
    return events


def events_to_phasic(events: list[dict], n: int, fs: float, tau0: float, tau1: float) -> np.ndarray:
    """Sum of Bateman responses a_k h(t - t_k) at the detected onsets (signal units)."""
    from scipy.signal import fftconvolve

    from .decompose import bateman

    imp = np.zeros(n)
    for e in events:
        imp[e["onset"]] += e["area"]
    h = bateman(np.arange(0, 10 * tau0, 1 / fs), tau0, tau1)
    return fftconvolve(imp, h)[:n]


def response_rates(rec: Recording, events: list[dict], clean: np.ndarray) -> dict:
    """Maximum-likelihood Poisson rate per stage type: K / T over clean windows."""
    st = window_stage(rec)
    out = {}
    for g in ("rest", "load"):
        sel = clean & np.array([stage_type(s) == g for s in st])
        T = sel.sum() * WIN / rec.fs
        win_of = np.array([e["onset"] // WIN for e in events], dtype=int)
        K = int(np.isin(win_of, np.where(sel)[0]).sum())
        out[g] = {"K": K, "T_s": float(T), "rate_per_s": (K / T) if T > 0 else float("nan")}
    return out


def fit_tonic(t: np.ndarray, fs: float, period_grid: np.ndarray) -> dict:
    """Least-squares fit of s0 + b*tau + A_s sin(2 pi tau / T + phi).

    The form is a project decision (slow trend plus one slow oscillation); a
    published source for the sinusoidal term has not been verified.
    """
    tau = np.arange(t.size) / fs
    best = None
    for T in period_grid:
        X = np.column_stack([np.ones_like(tau), tau, np.sin(2 * np.pi * tau / T), np.cos(2 * np.pi * tau / T)])
        coef, *_ = np.linalg.lstsq(X, t, rcond=None)
        resid = t - X @ coef
        sse = float(resid @ resid)
        if best is None or sse < best[0]:
            best = (sse, T, coef)
    sse, T, coef = best
    sst = float(((t - t.mean()) ** 2).sum())
    return {
        "s0": float(coef[0]),
        "b": float(coef[1]),
        "A_s": float(np.hypot(coef[2], coef[3])),
        "T_t": float(T),
        "r2": 1 - sse / sst if sst > 0 else float("nan"),
    }


def tonic_by_stage(rec: Recording, t: np.ndarray, cfg: dict) -> list[dict]:
    c = cfg["estimation"]["tonic"]
    grid = np.arange(c["period_min_s"], c["period_max_s"] + 1e-9, c["period_step_s"])
    out = []
    names = rec.stage
    starts = np.r_[0, np.where(names[1:] != names[:-1])[0] + 1]
    ends = np.r_[starts[1:], rec.n]
    for s, e in zip(starts, ends, strict=True):
        if (e - s) / rec.fs < 2 * c["period_max_s"]:
            continue  # too short to identify a 45-90 s oscillation
        fit = fit_tonic(t[s:e], rec.fs, grid)
        fit.update(stage=str(names[s]), stage_type=stage_type(str(names[s])))
        out.append(fit)
    return out


# ---------------------------------------------------------------- artifacts
def artifact_episodes(rec: Recording, sigma: float, cfg: dict) -> list[dict]:
    """Measure every run of majority-artifact windows and assign a class.

    Local baseline: straight line between the medians of the 2 s before and after
    the episode. Classification rules are fixed in configs (estimation.episodes)
    before the data are inspected.
    """
    c = cfg["estimation"]["episodes"]
    maj = rec.majority()
    y = rec.eda
    fs = rec.fs
    pad = int(round(c["baseline_s"] * fs))
    eps = []
    k = 0
    while k < maj.size:
        if not maj[k]:
            k += 1
            continue
        s_w = k
        while k < maj.size and maj[k]:
            k += 1
        s, e = s_w * WIN, k * WIN
        pre = y[max(0, s - pad) : s]
        post = y[e : min(y.size, e + pad)]
        if pre.size == 0 or post.size == 0:
            continue
        b0, b1 = float(np.median(pre)), float(np.median(post))
        base = np.linspace(b0, b1, e - s)
        dev = y[s:e] - base
        i_max = int(np.argmax(np.abs(dev)))
        A = float(dev[i_max])
        absd = np.abs(dev)
        above10 = np.where(absd >= 0.1 * abs(A))[0]
        above90 = np.where(absd >= 0.9 * abs(A))[0]
        rise = (above90[0] - above10[0]) / fs if above10.size and above90.size else 0.0
        ret = (above10[-1] - above90[-1]) / fs if above10.size and above90.size else 0.0
        seg = y[s:e] - y[s:e].mean()
        spec = np.abs(np.fft.rfft(seg)) ** 2
        f = np.fft.rfftfreq(seg.size, 1 / fs)
        pn = float(spec[f > 0.4].sum() / spec.sum()) if spec.sum() > 0 else 0.0
        shift = b1 - b0
        ep = {
            "pid": rec.pid,
            "start": s,
            "end": e,
            "duration_s": (e - s) / fs,
            "amp": A,
            "max_slope": float(np.max(np.abs(np.diff(y[s:e]))) * fs) if e - s > 1 else 0.0,
            "rise_s": float(rise),
            "return_s": float(ret),
            "pn": pn,
            "shift": float(shift),
            "stage": str(rec.stage[s]) if rec.stage is not None else "",
        }
        ep["cls"] = classify_episode(ep, y[s:e], sigma, c)
        eps.append(ep)
    return eps


def classify_episode(ep: dict, seg: np.ndarray, sigma: float, c: dict) -> str:
    """Rules fixed in advance (configs: estimation.episodes); order matters."""
    flat = np.diff(seg) == 0
    run = longest_run(flat) / 8.0
    if (seg.min() <= c["dropout_level"] or seg.max() >= c["saturation_level"]) or (
        run >= c["flat_min_s"] and (seg.min() <= c["dropout_level"] + 1e-9)
    ):
        return "C5"
    if abs(ep["amp"]) < c["min_amp"]:
        return "unclassified"
    if abs(ep["shift"]) > c["shift_sigma"] * sigma and abs(ep["shift"]) >= c["min_amp"]:
        return "C3"
    if ep["pn"] >= c["pn_walk"] and ep["stage"] in c["walk_stages"]:
        return "C4"
    return "C2" if ep["amp"] < 0 else "C1"


def longest_run(mask: np.ndarray) -> int:
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def template_parameters(ep: dict) -> dict:
    """Convert episode measurements into the parameters of the generator templates.

    For linear edges the 10-90 % time is 0.8 of the full edge time.
    """
    tr = max(ep["rise_s"] / 0.8, 1 / 64)
    tf = max(ep["return_s"] / 0.8, 1 / 64)
    cls = ep["cls"]
    if cls == "C1":
        return {"A": abs(ep["amp"]), "tau_r": tr, "tau_f": tf}
    if cls == "C2":
        w = tr
        return {"depth": abs(ep["amp"]), "w": w, "D": max(ep["duration_s"] - 2 * w, 0.0)}
    if cls == "C3":
        return {"delta": ep["shift"], "w": tr}
    if cls == "C5":
        return {"D": ep["duration_s"], "kind": "zero" if ep["amp"] < 0 else "saturation"}
    return {}


# ---------------------------------------------------------------- motion
def motion_events(rec: Recording, cfg: dict) -> tuple[list[dict], float]:
    """Events where the jerk exceeds the participant's own rest percentile.

    Runs of j above theta lasting >= min_len samples are kept first; kept runs
    closer than ``merge_s`` are then merged (merging isolated single-sample
    exceedances first would create spurious events).
    Returns the events (onset, end, intensity = max jerk) and theta.
    """
    c = cfg["estimation"]["motion"]
    j = jerk(rec.acc, rec.fs)
    rest = np.array([stage_type(str(s)) == "rest" for s in rec.stage])
    theta = float(np.percentile(j[rest] if rest.any() else j, c["percentile"]))
    above = j > theta
    runs = []
    k = 0
    while k < above.size:
        if above[k]:
            s = k
            while k < above.size and above[k]:
                k += 1
            if k - s >= c["min_len"]:  # length rule applies to raw runs, before merging
                runs.append([s, k])
        else:
            k += 1
    merged = []
    gap = int(round(c["merge_s"] * rec.fs))
    for r in runs:
        if merged and r[0] - merged[-1][1] < gap:
            merged[-1][1] = r[1]
        else:
            merged.append(r)
    ev = [{"onset": s, "end": e, "intensity": float(j[s:e].max())} for s, e in merged]
    return ev, theta


def utd_reference_point(episodes: list[dict], events: list[dict], fs: float, cfg: dict) -> dict:
    """Describe how the UTD labels relate to movement (not reality itself)."""
    lo, hi = cfg["estimation"]["reference"]["window_s"]
    on_e = np.array([e["onset"] for e in events]) / fs
    delays, linked = [], 0
    for ep in episodes:
        t = ep["start"] / fs
        d = t - on_e
        ok = d[(d >= lo) & (d <= hi)]
        if ok.size:
            linked += 1
            delays.append(float(ok.min()))
    followed = 0
    st = np.array([ep["start"] for ep in episodes]) / fs
    for t in on_e:
        d = st - t
        if np.any((d >= lo) & (d <= hi)):
            followed += 1
    return {
        "n_episodes": len(episodes),
        "n_motion_events": len(events),
        "share_episodes_after_motion": linked / len(episodes) if episodes else float("nan"),
        "share_motion_followed_by_artifact": followed / len(events) if events else float("nan"),
        "delays_s": delays,
    }


def walking_gain(rec: Recording, cfg: dict) -> dict:
    """Regress the fast part of EDA on the fast part of |a| in walking windows (class C4)."""
    c = cfg["estimation"]["c4"]
    k = int(round(c["smooth_s"] * rec.fs))
    kern = np.ones(k) / k
    mag = np.sqrt((rec.acc**2).sum(1))
    fe = rec.eda - np.convolve(rec.eda, kern, "same")
    fa_ = mag - np.convolve(mag, kern, "same")
    sel = np.array([str(s) in c["stages"] for s in rec.stage])
    if sel.sum() < 10 * k:
        return {"G": float("nan"), "r2": float("nan"), "n": int(sel.sum())}
    x, y = fa_[sel], fe[sel]
    G = float(np.dot(x, y) / np.dot(x, x))
    r2 = float(1 - ((y - G * x) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return {"G": G, "r2": r2, "n": int(sel.sum())}


# ---------------------------------------------------------------- plausibility
def plausibility_warnings(rates: dict, dec: dict, cfg: dict) -> list[str]:
    """Flags for values outside commonly reported ranges. Warnings, not gates.

    A response rate above ``max_rate_per_min`` (a project decision, set well
    above the few responses per minute usually described at rest; the exact
    literature range is still to be verified) suggests that the driver is
    absorbing artifacts or noise. tau0 at the edge of its grid means the optimum
    may lie outside it.
    """
    c = cfg["estimation"].get("plausibility", {"max_rate_per_min": 10.0})
    w = []
    for g, v in rates.items():
        if np.isfinite(v["rate_per_s"]) and 60 * v["rate_per_s"] > c["max_rate_per_min"]:
            w.append(f"{g} response rate {60 * v['rate_per_s']:.1f}/min above {c['max_rate_per_min']}/min")
    if dec["tau0_at_bound"]:
        w.append(f"tau0 = {dec['tau0']} at the edge of its grid")
    if dec.get("kkt") is not None and dec["kkt"] > 1e-3:
        w.append(f"decomposition not converged (KKT violation {dec['kkt']:.2e})")
    return w


# ---------------------------------------------------------------- one participant
def estimate_participant(rec: Recording, cfg: dict) -> dict:
    """Run steps 3.1-3.7 for one labelled recording."""
    clean = clean_windows(rec)
    sigma, sigma_source = noise_sigma(rec, clean, with_source=True)
    if not np.isfinite(sigma):
        raise ValueError(f"{rec.pid}: no clean windows to estimate the noise level")
    q = quantisation_step(rec.eda)
    q_acc = quantisation_step(rec.acc[:, 0])
    dec = decompose_participant(rec, clean, cfg)
    c = cfg["estimation"]
    ev = scr_events(dec["r"], rec.fs, dec["tau0"], c["cvxeda"]["tau1"], c["min_response_amp"])
    rates = response_rates(rec, ev, clean)
    # The generator rebuilds the signal as tonic + detected responses. cvxEDA's
    # phasic part also holds activity below the response threshold, so its tonic
    # alone sits too low (on AWW by several uS). The level the generator needs is
    # therefore the signal without the detected responses: t + r - r_events.
    r_ev = events_to_phasic(ev, rec.n, rec.fs, dec["tau0"], c["cvxeda"]["tau1"])
    tonic = tonic_by_stage(rec, dec["t"] + dec["r"] - r_ev, cfg)
    # noise floor for episode classification: the measured sigma can be 0 when the
    # noise is smaller than the quantisation step; rounding alone has sd q / sqrt(12)
    eps = artifact_episodes(rec, max(sigma, q / np.sqrt(12.0)), cfg)
    mev, theta = motion_events(rec, cfg)
    ref = utd_reference_point(eps, mev, rec.fs, cfg)
    return {
        "pid": rec.pid,
        "n_samples": rec.n,
        "n_clean_windows": int(clean.sum()),
        "sigma": sigma,
        "sigma_source": sigma_source,
        "q": q,
        "q_acc": q_acc,
        "tau0": dec["tau0"],
        "tau0_at_bound": dec["tau0_at_bound"],
        "tau0_scores": dec["tau0_scores"],
        "decomposition_backend": dec["backend"],
        "decomposition_kkt": dec["kkt"],
        "decomposition_iterations": dec["iterations"],
        "warnings": plausibility_warnings(rates, dec, cfg),
        "areas": [e["area"] for e in ev],
        "amps": [e["amp"] for e in ev],
        "rates": rates,
        "tonic": tonic,
        "episodes": eps,
        "motion_events": mev,
        "motion_theta": theta,
        "reference": ref,
        "c4": walking_gain(rec, cfg),
    }
