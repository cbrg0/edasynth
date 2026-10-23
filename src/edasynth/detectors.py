"""Untrained EDA-only detectors used as reference points (method step 7.3, goal 2).

None of them sees the accelerometer. They show which artifact classes are caught
without movement information; they are not part of the map.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import firwin, welch

from .features import WIN, windows


def _dilate(mask: np.ndarray, k: int) -> np.ndarray:
    if k <= 0:
        return mask
    return np.convolve(mask.astype(float), np.ones(2 * k + 1), "same") > 0


def to_windows(sample_mask: np.ndarray) -> np.ndarray:
    return windows(sample_mask).any(1).astype(int)


def slope_score(y: np.ndarray, fs: float) -> np.ndarray:
    """Continuous score: maximum |dy/dt| (forward difference) in each window."""
    v = np.zeros_like(y)
    v[1:] = np.abs(np.diff(y)) * fs
    return windows(v).max(1)


def kleckner(y: np.ndarray, fs: float, cfg: dict) -> np.ndarray:
    """Rules 1, 2 and 4 of Kleckner et al. (2018); rule 3 (temperature) omitted.

    The slope is computed on a zero-phase low-pass FIR filtered signal
    (cut-off 0.35 Hz, 64 s long, matching the duration of the authors' filter).
    """
    c = cfg["detectors"]["kleckner"]
    taps = int(round(c["fir_seconds"] * fs)) | 1
    h = firwin(taps, c["cutoff_hz"], fs=fs)
    yf = np.convolve(np.pad(y, taps // 2, mode="edge"), h, "valid")
    slope = np.zeros_like(yf)
    slope[1:] = np.diff(yf) * fs
    bad = (yf < c["level_min"]) | (yf > c["level_max"]) | (np.abs(slope) > c["slope_max"])
    bad = _dilate(bad, int(round(c["transition_s"] * fs)))
    return to_windows(bad)


def kocielnik(y: np.ndarray, fs: float, cfg: dict) -> np.ndarray:
    """Relative change faster than +20 %/s or -10 %/s (Kocielnik et al., cited by Taylor et al., 2015)."""
    c = cfg["detectors"]["kocielnik"]
    rel = np.zeros_like(y)
    base = np.where(np.abs(y[:-1]) > 1e-9, y[:-1], np.nan)
    rel[1:] = (np.diff(y) / base) * fs
    bad = (rel > c["rise_per_s"]) | (rel < -c["fall_per_s"])
    return to_windows(np.nan_to_num(bad, nan=True).astype(bool))


def pn_score(y: np.ndarray, fs: float) -> np.ndarray:
    """Share of power above 0.4 Hz in each window (Hossain et al., 2022, eq. 1)."""
    out = []
    for w in windows(y):
        f, pxx = welch(w - w.mean(), fs=fs, nperseg=WIN // 2, noverlap=WIN // 4)
        tot = pxx.sum()
        out.append(pxx[f > 0.4].sum() / tot if tot > 0 else 0.0)
    return np.asarray(out)


# ------------------------------------------------------------- Chen et al. (2015)
def swt_haar(x: np.ndarray, levels: int) -> list[np.ndarray]:
    """Undecimated (stationary) Haar transform, periodic extension; details d1..dJ."""
    a = np.asarray(x, float)
    out = []
    for j in range(levels):
        s = 2**j
        shifted = np.roll(a, -s)
        out.append((a - shifted) / np.sqrt(2))
        a = (a + shifted) / np.sqrt(2)
    return out


def _em_scale_mixture(d: np.ndarray, iters: int = 200):
    """EM for gamma N(0, s2) + (1 - gamma) N(0, c2 s2) with c2 > 1."""
    v = d**2
    s2 = np.median(v) / 0.4549 + 1e-12
    c2, g = 10.0, 0.9
    for _ in range(iters):
        p1 = g * np.exp(-v / (2 * s2)) / np.sqrt(s2)
        p2 = (1 - g) * np.exp(-v / (2 * c2 * s2)) / np.sqrt(c2 * s2)
        w = p1 / (p1 + p2 + 1e-300)
        g = np.clip(w.mean(), 1e-3, 1 - 1e-3)
        s2 = max((w * v).sum() / max(w.sum(), 1e-12), 1e-12)
        c2 = max(((1 - w) * v).sum() / max((1 - w).sum(), 1e-12) / s2, 1.0 + 1e-6)
    return g, s2, c2


def _mixture_quantile(prob: float, g: float, s2: float, c2: float) -> float:
    from scipy.stats import norm

    def cdf(x):
        return g * norm.cdf(x / np.sqrt(s2)) + (1 - g) * norm.cdf(x / np.sqrt(c2 * s2))

    lo, hi = -50 * np.sqrt(c2 * s2), 50 * np.sqrt(c2 * s2)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < prob:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def chen_swt(y: np.ndarray, fs: float, cfg: dict) -> np.ndarray:
    """Windows whose SWT coefficients fall outside the Gaussian-mixture thresholds.

    Haar, 8 levels, delta = 0.01, local windows of L = 400 s (Chen et al., 2015).
    Chen et al. remove such coefficients; marking a 5-s window as an artifact when
    at least one coefficient inside it was zeroed is our rule.
    """
    c = cfg["detectors"]["chen"]
    J, delta, L = c["levels"], c["delta"], int(c["window_s"] * fs)
    n = y.size
    m = 2**J
    npad = int(np.ceil(n / m) * m)
    x = np.pad(y, (0, npad - n), mode="symmetric")
    det = swt_haar(x, J)
    bad = np.zeros(npad, bool)
    for d in det:
        for s in range(0, npad, L):
            seg = d[s : s + L]
            if seg.size < 16:
                continue
            g, s2, c2 = _em_scale_mixture(seg)
            lo = _mixture_quantile(delta / 2, g, s2, c2)
            hi = _mixture_quantile(1 - delta / 2, g, s2, c2)
            bad[s : s + L] |= (seg < lo) | (seg > hi)
    return to_windows(bad[:n])


def run_all(y: np.ndarray, fs: float, cfg: dict) -> dict:
    """Scores (continuous) and flags (binary) of every detector for one recording."""
    return {
        "slope": ("score", slope_score(y, fs)),
        "pn": ("score", pn_score(y, fs)),
        "kleckner": ("flag", kleckner(y, fs, cfg)),
        "kocielnik": ("flag", kocielnik(y, fs, cfg)),
        "chen_swt": ("flag", chen_swt(y, fs, cfg)),
    }
