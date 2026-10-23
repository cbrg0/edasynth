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

WIN = 40

F0_NAMES = [
    f"eda_{sig}_{st}" for sig in ("amp", "d1", "d2", "w1hz", "w2hz", "w4hz") for st in ("mean", "sd", "max", "min")
]
_ACC_BLOCKS = [f"{ch}_{sig}" for ch in ("mag", "x", "y", "z") for sig in ("amp", "d1", "d2")] + [
    f"{ch}_{w}" for ch in ("mag", "x", "y", "z") for w in ("w1hz", "w2hz", "w4hz")
]
FA_NAMES = [f"acc_{b}_{st}" for b in _ACC_BLOCKS for st in ("mean", "sd", "max", "min")]
FA4_NAMES = FA_NAMES[:4]
FB_NAMES = ["B1_lagcorr", "B2_lag_s", "B3_cooccur", "B4_coherence"]


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


def haar_details(w: np.ndarray, levels: int = 3) -> list[np.ndarray]:
    """Haar DWT detail coefficients [cD1, cD2, cD3] along the last axis.

    Uses the PyWavelets sign convention cD[k] = (x[2k] - x[2k+1]) / sqrt(2).
    For a 40-sample window every level has even length (40, 20, 10), so no
    boundary extension is involved and the result equals
    ``pywt.wavedec(w, 'haar', level=3)`` exactly (checked in tests when pywt is
    installed).
    """
    out = []
    a = np.asarray(w, dtype=float)
    for _ in range(levels):
        if a.shape[-1] % 2:
            raise ValueError("odd length: boundary handling would be needed")
        even, odd = a[..., 0::2], a[..., 1::2]
        out.append((even - odd) / np.sqrt(2.0))
        a = (even + odd) / np.sqrt(2.0)
    return out


def windows(x: np.ndarray, win: int = WIN) -> np.ndarray:
    """Reshape a 1-D signal (or (n, k) array) into (n_windows, win[, k])."""
    nw = x.shape[0] // win
    return x[: nw * win].reshape(nw, win, *x.shape[1:])


def _block24(w: np.ndarray, deriv: str) -> np.ndarray:
    d1, d2 = derivatives(w, deriv)
    c1, c2, c3 = haar_details(w)
    return np.concatenate([stats4(w), stats4(d1), stats4(d2), stats4(c3), stats4(c2), stats4(c1)], -1)


# ---------------------------------------------------------------- feature sets
def f0(eda: np.ndarray, deriv: str = "central") -> np.ndarray:
    """24 EDA features per window, shape (n_windows, 24)."""
    return _block24(windows(eda), deriv)


def fa(acc: np.ndarray, deriv: str = "central") -> np.ndarray:
    """96 accelerometer features per window in the order of Zhang et al."""
    mag = np.sqrt((acc**2).sum(1))
    chans = [windows(mag), windows(acc[:, 0]), windows(acc[:, 1]), windows(acc[:, 2])]
    stat_part, wave_part = [], []
    for w in chans:
        d1, d2 = derivatives(w, deriv)
        stat_part += [stats4(w), stats4(d1), stats4(d2)]
        c1, c2, c3 = haar_details(w)
        wave_part += [stats4(c3), stats4(c2), stats4(c1)]
    return np.concatenate(stat_part + wave_part, -1)


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
    q_eda: float,
    q_acc: float,
    max_lag_s: float = 5.0,
    context_s: float = 5.0,
    median_s: float = 60.0,
    coh_band: tuple[float, float] = (0.5, 3.0),
    coh_seg_s: float = 2.5,
) -> np.ndarray:
    """Four causal coupling features per window, shape (n_windows, 4).

    B1  max over lags 0..L of corr(j[n - lag], v[n]) for n in the window,
        where v[n] = fs * |y[n] - y[n-1]| (speed of EDA change) and j is the jerk.
    B2  the lag of B1, in seconds.
    B3  [max_w v / (median_60s(v) + q_eda*fs)] * [max_C j / (median_60s(j) + q_acc*fs)];
        medians over the 60 s preceding the window (the window itself if none).
    B4  mean magnitude-squared coherence of |a| and y in ``coh_band`` over the
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
    lag = int(round(max_lag_s * fs))
    ctx = int(round(context_s * fs))
    med = int(round(median_s * fs))
    seg = int(round(coh_seg_s * fs))
    out = np.zeros((nw, 4))
    for i in range(nw):
        s, e = i * WIN, (i + 1) * WIN
        vw = v[s:e]
        # B1, B2: lagged correlation using only past or current jerk
        best, best_lag = 0.0, 0
        if vw.std() > 0:
            for ell in range(0, min(lag, s) + 1):
                jw = j[s - ell : e - ell]
                if jw.std() > 0:
                    r = float(np.corrcoef(jw, vw)[0, 1])
                    if r > best:
                        best, best_lag = r, ell
        out[i, 0] = best
        out[i, 1] = best_lag / fs
        # B3
        ps = max(0, s - med)
        past_v = v[ps:s] if s > ps else vw
        cs = max(0, s - ctx)
        past_j = j[ps:s] if s > ps else j[cs:e]
        out[i, 2] = (vw.max() / (np.median(past_v) + q_eda * fs)) * (j[cs:e].max() / (np.median(past_j) + q_acc * fs))
        # B4
        xa, ya = mag[cs:e], eda[cs:e]
        if xa.std() > 0 and ya.std() > 0:
            with np.errstate(invalid="ignore", divide="ignore"):
                f, cxy = coherence(xa, ya, fs=fs, nperseg=seg, noverlap=seg // 2)
            band = (f >= coh_band[0]) & (f <= coh_band[1])
            vals = cxy[band][np.isfinite(cxy[band])]
            # a constant segment has no spectrum: coherence undefined -> 0 (no shared change)
            out[i, 3] = float(vals.mean()) if vals.size else 0.0
    return out


def feature_matrix(
    eda: np.ndarray,
    acc: np.ndarray,
    fs: float,
    q_eda: float,
    q_acc: float,
    deriv: str = "central",
    fb_kwargs: dict | None = None,
) -> dict[str, np.ndarray]:
    """All feature sets for one recording."""
    return {
        "F0": f0(eda, deriv),
        "FA": fa(acc, deriv),
        "FA4": fa4(acc),
        "FB": fb(eda, acc, fs, q_eda, q_acc, **(fb_kwargs or {})),
    }
