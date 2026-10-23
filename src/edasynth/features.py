"""Window features.

All features use non-overlapping 5-s windows (40 samples at 8 Hz), the unit used
by Taylor et al. (2015) and Zhang et al. (2017). The same functions are applied to
real and synthetic recordings.

Feature sets
------------
F0  : 24 EDA features of Zhang et al. (2017): mean, SD, max, min of the signal,
      its first and second derivatives and Haar detail coefficients at levels 3, 2, 1
      (1, 2 and 4 Hz), in that order.
FA  : 96 accelerometer features of Zhang et al. (2017), same statistics for the
      acceleration magnitude and each axis (order as in their "Feature Details.txt").
FA4 : the first 4 of FA (mean, SD, max, min of acceleration magnitude).
FB  : 4 coupling features proposed here (B1-B4), all causal: window w uses only the
      5 s before it and w itself.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import coherence

from .config import key

WIN = 40

f0Names = [
    key("eda", sig, st) for sig in ("amp", "d1", "d2", "w1hz", "w2hz", "w4hz") for st in ("mean", "sd", "max", "min")
]
accBlocks = [key(ch, sig) for ch in ("mag", "x", "y", "z") for sig in ("amp", "d1", "d2")] + [
    key(ch, w) for ch in ("mag", "x", "y", "z") for w in ("w1hz", "w2hz", "w4hz")
]
faNames = [key("acc", b, st) for b in accBlocks for st in ("mean", "sd", "max", "min")]
fa4Names = faNames[:4]
fbNames = ["B1Lagcorr", "B2LagS", "B3Cooccur", "B4Coherence"]


# ---------------------------------------------------------------- primitives
def stats4(x: np.ndarray) -> np.ndarray:
    """Mean, population SD, max, min along the last axis."""
    return np.stack([x.mean(-1), x.std(-1), x.max(-1), x.min(-1)], axis=-1)


def derivatives(w: np.ndarray, mode: str = "central") -> tuple[np.ndarray, np.ndarray]:
    """First and second differences along the last axis.

    ``central``: (x[n+1] - x[n-1]) / 2, as in the published script of Zhang et al.
    ``forward``: x[n+1] - x[n], which reproduces their stored feature matrices.
    """
    if mode == "central":
        d1 = (w[..., 2:] - w[..., :-2]) / 2.0
    elif mode == "forward":
        d1 = np.diff(w, axis=-1)
    else:
        raise ValueError(mode)
    d2 = w[..., 2:] - 2 * w[..., 1:-1] + w[..., :-2]
    return d1, d2


def haarDetails(w: np.ndarray, levels: int = 3) -> list[np.ndarray]:
    """Haar DWT detail coefficients [cD1, cD2, cD3] along the last axis.

    Uses the PyWavelets sign convention cD[k] = (x[2k] - x[2k+1]) / sqrt(2).
    For a 40-sample window every level has even length (40, 20, 10), so no
    boundary extension is involved and the result equals
    ``pywt.wavedec(w, 'haar', level=3)`` exactly (checked in tests when pywt is
    installed).
    """
    out = []
    a = np.asarray(w, dtype=float)
    for iteration in range(levels):
        if a.shape[-1] % 2:
            raise ValueError("odd length: boundary handling would be needed")
        even, odd = a[..., 0::2], a[..., 1::2]
        out.append((even - odd) / np.sqrt(2.0))
        a = (even + odd) / np.sqrt(2.0)
    return out


def windows(x: np.ndarray, win: int = WIN) -> np.ndarray:
    """Reshape a 1-D signal (or (n, k) array) into (nWindows, win[, k])."""
    nw = x.shape[0] // win
    return x[: nw * win].reshape(nw, win, *x.shape[1:])


def block24(w: np.ndarray, deriv: str) -> np.ndarray:
    d1, d2 = derivatives(w, deriv)
    c1, c2, c3 = haarDetails(w)
    return np.concatenate([stats4(w), stats4(d1), stats4(d2), stats4(c3), stats4(c2), stats4(c1)], -1)


# ---------------------------------------------------------------- feature sets
def f0(eda: np.ndarray, deriv: str = "central") -> np.ndarray:
    """24 EDA features per window, shape (nWindows, 24)."""
    return block24(windows(eda), deriv)


def fa(acc: np.ndarray, deriv: str = "central") -> np.ndarray:
    """96 accelerometer features per window in the order of Zhang et al."""
    mag = np.sqrt((acc**2).sum(1))
    chans = [windows(mag), windows(acc[:, 0]), windows(acc[:, 1]), windows(acc[:, 2])]
    statPart, wavePart = [], []
    for w in chans:
        d1, d2 = derivatives(w, deriv)
        statPart += [stats4(w), stats4(d1), stats4(d2)]
        c1, c2, c3 = haarDetails(w)
        wavePart += [stats4(c3), stats4(c2), stats4(c1)]
    return np.concatenate(statPart + wavePart, -1)


def fa4(acc: np.ndarray) -> np.ndarray:
    """Mean, SD, max, min of the acceleration magnitude per window."""
    return stats4(windows(np.sqrt((acc**2).sum(1))))


def jerk(acc: np.ndarray, fs: float) -> np.ndarray:
    """j[n] = fs * ||a[n] - a[n-1]|| in g/s, with j[0] = 0."""
    j = np.zeros(acc.shape[0])
    j[1:] = fs * np.sqrt((np.diff(acc, axis=0) ** 2).sum(1))
    return j


def fb(
    eda: np.ndarray,
    acc: np.ndarray,
    fs: float,
    qEda: float,
    qAcc: float,
    maxLagS: float = 5.0,
    contextS: float = 5.0,
    medianS: float = 60.0,
    cohBand: tuple[float, float] = (0.5, 3.0),
    cohSegS: float = 2.5,
) -> np.ndarray:
    """Four causal coupling features per window, shape (nWindows, 4).

    B1  max over lags 0..L of corr(j[n - lag], v[n]) for n in the window,
        where v[n] = fs * |y[n] - y[n-1]| (speed of EDA change) and j is the jerk.
    B2  the lag of B1, in seconds.
    B3  [max v in the window / (median v over 60 s + qEda*fs)] *
        [max j in the context / (median j over 60 s + qAcc*fs)];
        medians over the 60 s preceding the window (the window itself if none).
    B4  mean magnitude-squared coherence of |a| and y in ``cohBand`` over the
        context C(w) = 5 s before the window plus the window (Welch, 2.5-s segments,
        50 % overlap).
    The quantisation terms keep the ratios finite on flat quantised signals;
    B4 is set to 0 where a signal is constant and coherence is undefined.
    """
    n = eda.size
    nw = n // WIN
    v = np.zeros(n)
    v[1:] = fs * np.abs(np.diff(eda))
    j = jerk(acc, fs)
    mag = np.sqrt((acc**2).sum(1))
    lag = int(round(maxLagS * fs))
    ctx = int(round(contextS * fs))
    med = int(round(medianS * fs))
    seg = int(round(cohSegS * fs))
    out = np.zeros((nw, 4))
    for i in range(nw):
        s, e = i * WIN, (i + 1) * WIN
        vw = v[s:e]
        # B1, B2: lagged correlation using only past or current jerk
        best, bestLag = 0.0, 0
        if vw.std() > 0:
            for ell in range(0, min(lag, s) + 1):
                jw = j[s - ell : e - ell]
                if jw.std() > 0:
                    r = float(np.corrcoef(jw, vw)[0, 1])
                    if r > best:
                        best, bestLag = r, ell
        out[i, 0] = best
        out[i, 1] = bestLag / fs
        # B3
        ps = max(0, s - med)
        pastV = v[ps:s] if s > ps else vw
        cs = max(0, s - ctx)
        pastJ = j[ps:s] if s > ps else j[cs:e]
        out[i, 2] = (vw.max() / (np.median(pastV) + qEda * fs)) * (j[cs:e].max() / (np.median(pastJ) + qAcc * fs))
        # B4
        xa, ya = mag[cs:e], eda[cs:e]
        if xa.std() > 0 and ya.std() > 0:
            with np.errstate(invalid="ignore", divide="ignore"):
                f, cxy = coherence(xa, ya, fs=fs, nperseg=seg, noverlap=seg // 2)
            band = (f >= cohBand[0]) & (f <= cohBand[1])
            vals = cxy[band][np.isfinite(cxy[band])]
            # a constant segment has no spectrum: coherence undefined -> 0 (no shared change)
            out[i, 3] = float(vals.mean()) if vals.size else 0.0
    return out


def featureMatrix(
    eda: np.ndarray,
    acc: np.ndarray,
    fs: float,
    qEda: float,
    qAcc: float,
    deriv: str = "central",
    fbKwargs: dict | None = None,
) -> dict[str, np.ndarray]:
    """All feature sets for one recording."""
    return {
        "F0": f0(eda, deriv),
        "FA": fa(acc, deriv),
        "FA4": fa4(acc),
        "FB": fb(eda, acc, fs, qEda, qAcc, **(fbKwargs or {})),
    }
