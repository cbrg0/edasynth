"""This file is the third step of the study, measuring the real recordings. It learns from
the UTD recordings every property that the generator later reproduces, so that as little
as possible in the synthetic data is invented.

Each measurement follows a source wherever one exists. A response counts only when it is
at least 0.05 microsiemens high, the threshold accepted in the literature (Kleckner et
al., 2018; Boucsein et al., 2012, as summarised by Posada-Quintero & Chon, 2020). Clean
windows are taken only with a clean neighbour on each side, after the margin around
invalid data of Kleckner et al. (2018). The share of fast oscillation in an artifact is
measured above 0.4 Hz, as Hossain et al. (2022) did. The usual estimate of noise (Donoho
& Johnstone, 1994) gave zero on these recordings because the noise is smaller than the
device's smallest step, so the noise is instead taken directly from the part of the
signal that the decomposition leaves unexplained. Where no source exists, as for the
rules that sort artifact episodes into types or the threshold that defines a movement,
the values are decisions of this study, written into the settings before the main
results and listed in study-description/METHOD.md.

The outcome is, for every participant, the noise, the rates and sizes of responses, the
course of the slow level, the measured shape of each artifact and the moments of
movement.
"""

from __future__ import annotations

import numpy as np

from .data import WIN, Recording, stageType
from .decompose import decompose, hMax
from .features import haarDetails, jerk, windows


def cleanWindows(rec: Recording) -> np.ndarray:
    unanClean = rec.labels3.sum(1) == 0
    ok = unanClean.copy()
    ok[1:] &= unanClean[:-1]
    ok[:-1] &= unanClean[1:]
    return ok


def windowMaskToSamples(mask: np.ndarray, n: int) -> np.ndarray:
    m = np.zeros(n, bool)
    m[: mask.size * WIN] = np.repeat(mask, WIN)
    return m


def windowStage(rec: Recording) -> np.ndarray:
    return rec.stage[: rec.nWindows * WIN : WIN]


def noiseSigma(rec: Recording, clean: np.ndarray, withSource: bool = False):
    st = windowStage(rec)
    sel = clean & np.array([stageType(s) == "rest" for s in st])
    source = "cleanRest"
    if not sel.any():
        sel, source = clean.copy(), "cleanAllFallback"
    w = windows(rec.eda)[sel]
    if w.size == 0:
        sigma, source = float("nan"), "none"
    else:
        d1 = haarDetails(w)[0]
        sigma = float(np.median(np.abs(d1)) / 0.6745)
    return (sigma, source) if withSource else sigma


def residualBlocks(rec: Recording, residual: np.ndarray, clean: np.ndarray, blockS: float) -> dict:
    st = windowStage(rec)
    sel = clean & np.array([stageType(s) == "rest" for s in st])
    source = "cleanRest"
    if not sel.any():
        sel, source = clean.copy(), "cleanAllFallback"
    L = int(round(blockS * rec.fs))
    if L < 1 or WIN % L != 0:
        raise ValueError(f"noise block of {L} samples must divide the {WIN}-sample window")
    w = windows(residual)[sel]
    if w.size == 0:
        return {"blocks": np.zeros((0, L)), "source": "none", "sd": float("nan")}
    b = w.reshape(-1, L)
    b = b - b.mean(1, keepdims=True)
    return {"blocks": b, "source": source, "sd": float(b.std())}


def donorLevelCheck(eda: np.ndarray, cfg: dict) -> dict:
    c = cfg["estimation"]["donor"]
    share = float(np.mean(eda < c["minLevel"]))
    return {"shareBelowMinLevel": share, "validDonor": bool(share <= c["maxShareBelow"])}


def quantisationStep(x: np.ndarray, decimals: int = 9, minShare: float = 0.01) -> float:
    d = np.round(np.abs(np.diff(x)), decimals)
    d = d[d > 0]
    if d.size == 0:
        return 0.0
    vals, counts = np.unique(d, return_counts=True)
    ok = vals[counts >= minShare * d.size]
    return float(ok.min()) if ok.size else float(vals[np.argmax(counts)])


def decomposeParticipant(rec: Recording, clean: np.ndarray, cfg: dict) -> dict:
    c = cfg["estimation"]["cvxeda"]
    y = rec.eda.astype(float)
    mu, sd = float(y.mean()), float(y.std())
    z = (y - mu) / sd
    mask = windowMaskToSamples(clean, rec.n)
    best, pWarm = None, None
    grid = np.round(np.arange(c["tau0Min"], c["tau0Max"] + 1e-9, c["tau0Step"]), 3)
    scores = {}
    for tau0 in grid:
        out = decompose(
            z,
            1.0 / rec.fs,
            float(tau0),
            c["tau1"],
            c["alpha"],
            c["gamma"],
            c["knotS"],
            backend=c["backend"],
            p0=pWarm,
        )
        pWarm = out["p"]
        err = float(out["obj"])
        scores[float(tau0)] = {"objective": err, "residualClean": float(np.linalg.norm(out["e"][mask]))}
        if best is None or err < best[1]:
            best = (float(tau0), err, out)
    tau0, bestScore, out = best
    return {
        "tau0": tau0,
        "tau0Scores": scores,
        "kkt": out.get("kkt"),
        "iterations": out.get("iterations"),
        "r": sd * out["r"],
        "p": sd * out["p"],
        "t": sd * out["t"] + mu,
        "backend": out["backend"],
        "tau0AtBound": tau0 in (grid[0], grid[-1]),
    }


def scrEvents(r: np.ndarray, fs: float, tau0: float, tau1: float, minAmp: float) -> list[dict]:
    from scipy.signal import find_peaks

    pk, prop = find_peaks(r, prominence=minAmp)
    hm = hMax(tau0, tau1)
    back = int(round(5 * fs))
    events = []
    for k, A in zip(pk, prop["prominences"], strict=True):
        s0 = max(0, k - back)
        seg = r[s0 : k + 1]
        onset = s0 + seg.size - 1 - int(np.argmin(seg[::-1]))
        events.append({"onset": int(onset), "peak": int(k), "area": float(A / hm), "amp": float(A)})
    return events


def eventsToPhasic(events: list[dict], n: int, fs: float, tau0: float, tau1: float) -> np.ndarray:
    from scipy.signal import fftconvolve

    from .decompose import bateman

    imp = np.zeros(n)
    for e in events:
        imp[e["onset"]] += e["area"]
    h = bateman(np.arange(0, 10 * tau0, 1 / fs), tau0, tau1)
    return fftconvolve(imp, h)[:n]


def responseRates(rec: Recording, events: list[dict], clean: np.ndarray) -> dict:
    st = windowStage(rec)
    out = {}
    for g in ("rest", "load"):
        sel = clean & np.array([stageType(s) == g for s in st])
        T = sel.sum() * WIN / rec.fs
        winOf = np.array([e["onset"] // WIN for e in events], dtype=int)
        K = int(np.isin(winOf, np.where(sel)[0]).sum())
        out[g] = {"K": K, "TS": float(T), "ratePerS": (K / T) if T > 0 else float("nan")}
    return out


def fitTonic(t: np.ndarray, fs: float, periodGrid: np.ndarray) -> dict:
    tau = np.arange(t.size) / fs
    best = None
    for T in periodGrid:
        X = np.column_stack([np.ones_like(tau), tau, np.sin(2 * np.pi * tau / T), np.cos(2 * np.pi * tau / T)])
        coef, *rest = np.linalg.lstsq(X, t, rcond=None)
        resid = t - X @ coef
        sse = float(resid @ resid)
        if best is None or sse < best[0]:
            best = (sse, T, coef)
    sse, T, coef = best
    sst = float(((t - t.mean()) ** 2).sum())
    return {
        "s0": float(coef[0]),
        "b": float(coef[1]),
        "AS": float(np.hypot(coef[2], coef[3])),
        "TT": float(T),
        "r2": 1 - sse / sst if sst > 0 else float("nan"),
    }


def tonicByStage(rec: Recording, t: np.ndarray, cfg: dict) -> list[dict]:
    c = cfg["estimation"]["tonic"]
    grid = np.arange(c["periodMinS"], c["periodMaxS"] + 1e-9, c["periodStepS"])
    out = []
    names = rec.stage
    starts = np.r_[0, np.where(names[1:] != names[:-1])[0] + 1]
    ends = np.r_[starts[1:], rec.n]
    for s, e in zip(starts, ends, strict=True):
        if (e - s) / rec.fs < 2 * c["periodMaxS"]:
            continue
        fit = fitTonic(t[s:e], rec.fs, grid)
        fit.update(stage=str(names[s]), stageType=stageType(str(names[s])), durationS=(e - s) / rec.fs)
        out.append(fit)
    return out


def artifactEpisodes(rec: Recording, sigma: float, cfg: dict) -> list[dict]:
    c = cfg["estimation"]["episodes"]
    maj = rec.majority()
    y = rec.eda
    fs = rec.fs
    pad = int(round(c["baselineS"] * fs))
    eps = []
    k = 0
    while k < maj.size:
        if not maj[k]:
            k += 1
            continue
        sW = k
        while k < maj.size and maj[k]:
            k += 1
        s, e = sW * WIN, k * WIN
        pre = y[max(0, s - pad) : s]
        post = y[e : min(y.size, e + pad)]
        if pre.size == 0 or post.size == 0:
            continue
        b0, b1 = float(np.median(pre)), float(np.median(post))
        base = np.linspace(b0, b1, e - s)
        dev = y[s:e] - base
        iMax = int(np.argmax(np.abs(dev)))
        A = float(dev[iMax])
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
        lvl = max(b0, c["levelFloor"])
        ep = {
            "level": b0,
            "levelAfter": b1,
            "pid": rec.pid,
            "start": s,
            "end": e,
            "durationS": (e - s) / fs,
            "amp": A,
            "maxSlope": float(np.max(np.abs(np.diff(y[s:e]))) * fs) if e - s > 1 else 0.0,
            "riseS": float(rise),
            "returnS": float(ret),
            "pn": pn,
            "shift": float(shift),
            "stage": str(rec.stage[s]) if rec.stage is not None else "",
            "traceDev": [float(v) for v in dev / lvl],
            "traceStep": [float(v) for v in (y[s:e] - b0) / lvl],
        }
        ep["cls"] = classifyEpisode(ep, y[s:e], sigma, c)
        eps.append(ep)
    return eps


def classifyEpisode(ep: dict, seg: np.ndarray, sigma: float, c: dict) -> str:
    flat = np.diff(seg) == 0
    run = longestRun(flat) / 8.0
    if (seg.min() <= c["dropoutLevel"] or seg.max() >= c["saturationLevel"]) or (
        run >= c["flatMinS"] and (seg.min() <= c["dropoutLevel"] + 1e-9)
    ):
        return "C5"
    if abs(ep["amp"]) < c["minAmp"]:
        return "unclassified"
    if abs(ep["shift"]) > c["shiftSigma"] * sigma and abs(ep["shift"]) >= c["minAmp"]:
        return "C3"
    if ep["pn"] >= c["pnWalk"] and ep["stage"] in c["walkStages"]:
        return "C4"
    return "C2" if ep["amp"] < 0 else "C1"


def longestRun(mask: np.ndarray) -> int:
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def templateParameters(ep: dict, levelFloor: float) -> dict:
    tr = max(ep["riseS"] / 0.8, 1 / 64)
    tf = max(ep["returnS"] / 0.8, 1 / 64)
    cls = ep["cls"]
    lvl = max(ep["level"], levelFloor)
    if cls == "C1":
        return {"relA": abs(ep["amp"]) / lvl, "tauR": tr, "tauF": tf, "trace": ep.get("traceDev")}
    if cls == "C2":
        w = tr
        return {
            "relDepth": abs(ep["amp"]) / lvl,
            "w": w,
            "D": max(ep["durationS"] - 2 * w, 0.0),
            "trace": ep.get("traceDev"),
        }
    if cls == "C3":
        return {"relDelta": ep["shift"] / lvl, "w": tr, "trace": ep.get("traceStep")}
    if cls == "C5":
        return {"D": ep["durationS"], "kind": "zero" if ep["amp"] < 0 else "saturation"}
    return {}


def motionEvents(rec: Recording, cfg: dict) -> tuple[list[dict], float]:
    c = cfg["estimation"]["motion"]
    j = jerk(rec.acc, rec.fs)
    rest = np.array([stageType(str(s)) == "rest" for s in rec.stage])
    theta = float(np.percentile(j[rest] if rest.any() else j, c["percentile"]))
    above = j > theta
    runs = []
    k = 0
    while k < above.size:
        if above[k]:
            s = k
            while k < above.size and above[k]:
                k += 1
            if k - s >= c["minLen"]:
                runs.append([s, k])
        else:
            k += 1
    merged = []
    gap = int(round(c["mergeS"] * rec.fs))
    for r in runs:
        if merged and r[0] - merged[-1][1] < gap:
            merged[-1][1] = r[1]
        else:
            merged.append(r)
    ev = [{"onset": s, "end": e, "intensity": float(j[s:e].max())} for s, e in merged]
    return ev, theta


def utdReferencePoint(episodes: list[dict], events: list[dict], fs: float, cfg: dict) -> dict:
    lo, hi = cfg["estimation"]["reference"]["windowS"]
    onE = np.array([e["onset"] for e in events]) / fs
    delays, linked = [], 0
    for ep in episodes:
        t = ep["start"] / fs
        d = t - onE
        ok = d[(d >= lo) & (d <= hi)]
        if ok.size:
            linked += 1
            delays.append(float(ok.min()))
    followed = 0
    st = np.array([ep["start"] for ep in episodes]) / fs
    for t in onE:
        d = st - t
        if np.any((d >= lo) & (d <= hi)):
            followed += 1
    return {
        "nEpisodes": len(episodes),
        "nMotionEvents": len(events),
        "shareEpisodesAfterMotion": linked / len(episodes) if episodes else float("nan"),
        "shareMotionFollowedByArtifact": followed / len(events) if events else float("nan"),
        "delaysS": delays,
    }


def walkingGain(rec: Recording, cfg: dict) -> dict:
    c = cfg["estimation"]["c4"]
    k = int(round(c["smoothS"] * rec.fs))
    kern = np.ones(k) / k
    mag = np.sqrt((rec.acc**2).sum(1))
    fe = rec.eda - np.convolve(rec.eda, kern, "same")
    fastAcc = mag - np.convolve(mag, kern, "same")
    sel = np.array([str(s) in c["stages"] for s in rec.stage])
    if sel.sum() < 10 * k:
        return {"G": float("nan"), "r2": float("nan"), "n": int(sel.sum())}
    x, y = fastAcc[sel], fe[sel]
    G = float(np.dot(x, y) / np.dot(x, x))
    r2 = float(1 - ((y - G * x) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    return {"G": G, "r2": r2, "n": int(sel.sum())}


def cleanEvents(events: list[dict], clean: np.ndarray) -> list[dict]:
    last = clean.size - 1
    return [e for e in events if clean[min(e["onset"] // WIN, last)] and clean[min(e["peak"] // WIN, last)]]


def plausibilityWarnings(rates: dict, dec: dict, cfg: dict, amps=(), signalRange: float = float("inf")) -> list[str]:
    c = cfg["estimation"].get("plausibility", {"maxRatePerMin": 10.0})
    w = []
    for g, v in rates.items():
        if np.isfinite(v["ratePerS"]) and 60 * v["ratePerS"] > c["maxRatePerMin"]:
            w.append(f"{g} response rate {60 * v['ratePerS']:.1f}/min above {c['maxRatePerMin']}/min")
    if dec["tau0AtBound"]:
        w.append(f"tau0 = {dec['tau0']} at the edge of its grid")
    big = int(np.sum(np.asarray(amps, float) > signalRange))
    if big:
        w.append(f"{big} response height(s) larger than the range of the recording ({signalRange:.2f})")
    if dec.get("kkt") is not None and dec["kkt"] > 1e-3:
        w.append(f"decomposition not converged (KKT violation {dec['kkt']:.2e})")
    return w


def estimateParticipant(rec: Recording, cfg: dict) -> dict:
    clean = cleanWindows(rec)
    sigma, sigmaSource = noiseSigma(rec, clean, withSource=True)
    if not np.isfinite(sigma):
        raise ValueError(f"{rec.pid}: no clean windows to estimate the noise level")
    q = quantisationStep(rec.eda)
    qAcc = quantisationStep(rec.acc[:, 0])
    dec = decomposeParticipant(rec, clean, cfg)
    c = cfg["estimation"]
    ev = scrEvents(dec["r"], rec.fs, dec["tau0"], c["cvxeda"]["tau1"], c["minResponseAmp"])
    rates = responseRates(rec, ev, clean)
    sized = cleanEvents(ev, clean)
    if not sized:
        sized = [e for e in ev if clean[min(e["onset"] // WIN, clean.size - 1)]] or ev
    rEv = eventsToPhasic(ev, rec.n, rec.fs, dec["tau0"], c["cvxeda"]["tau1"])
    tonic = tonicByStage(rec, dec["t"] + dec["r"] - rEv, cfg)
    eps = artifactEpisodes(rec, max(sigma, q / np.sqrt(12.0)), cfg)
    resid = residualBlocks(rec, rec.eda - dec["t"] - dec["r"], clean, c["noise"]["blockS"])
    mev, theta = motionEvents(rec, cfg)
    ref = utdReferencePoint(eps, mev, rec.fs, cfg)
    return {
        "pid": rec.pid,
        "nSamples": rec.n,
        "nCleanWindows": int(clean.sum()),
        "sigma": sigma,
        "sigmaSource": sigmaSource,
        "residual": resid["blocks"],
        "residualSource": resid["source"],
        "residualSd": resid["sd"],
        **donorLevelCheck(rec.eda, cfg),
        "q": q,
        "qAcc": qAcc,
        "tau0": dec["tau0"],
        "tau0AtBound": dec["tau0AtBound"],
        "tau0Scores": dec["tau0Scores"],
        "decompositionBackend": dec["backend"],
        "decompositionKkt": dec["kkt"],
        "decompositionIterations": dec["iterations"],
        "warnings": plausibilityWarnings(rates, dec, cfg, [e["amp"] for e in sized], float(np.ptp(rec.eda))),
        "areas": [e["area"] for e in sized],
        "amps": [e["amp"] for e in sized],
        "nResponsesAll": len(ev),
        "rates": rates,
        "tonic": tonic,
        "episodes": eps,
        "motionEvents": mev,
        "motionTheta": theta,
        "reference": ref,
        "c4": walkingGain(rec, cfg),
    }
