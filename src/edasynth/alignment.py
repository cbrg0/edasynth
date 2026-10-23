"""Linking the expert labels (Zhang et al. repository) to the raw UTD signals.

Labels and signals were published separately and no document states which
window belongs to which piece of signal. Two checks are applied per participant:
1. the number of complete 5-s windows equals the number of labelled windows of
   group g = subject - 1;
2. the four amplitude features recomputed from the raw signal and z-scored within
   the participant (exactly as in Zhang's script) correlate >= 0.999 with the
   first four columns of the published matrix. If not, shifts of 1-10 windows and
   0-39 samples are searched; a participant that still fails is excluded.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .data import WIN
from .features import stats4, windows


def amplitude_z(eda: np.ndarray, offset: int = 0, n_windows: int | None = None) -> np.ndarray:
    w = windows(eda[offset:])
    if n_windows is not None:
        w = w[:n_windows]
    return stats.zscore(stats4(w), axis=0)


def check_participant(eda: np.ndarray, stored4: np.ndarray, threshold: float) -> dict:
    """Return window counts, best shift and per-feature correlations."""
    n_lab = stored4.shape[0]
    n_raw = eda.size // WIN

    def corr_at(off_samples: int, off_windows: int):
        z = amplitude_z(eda, off_samples + off_windows * WIN)
        m = min(z.shape[0], n_lab)
        if m < 10:
            return np.full(4, -1.0)
        z = stats.zscore(stats4(windows(eda[off_samples + off_windows * WIN :])[:m]), axis=0)
        ref = stats.zscore(stored4[:m], axis=0)
        return np.array([np.corrcoef(z[:, k], ref[:, k])[0, 1] for k in range(4)])

    r0 = corr_at(0, 0)
    best = {"shift_samples": 0, "shift_windows": 0, "r": r0}
    if r0.min() < threshold:
        for ow in range(-10, 11):
            for os_ in range(0, WIN):
                if ow * WIN + os_ < 0:
                    continue
                r = corr_at(os_, ow)
                if r.min() > best["r"].min():
                    best = {"shift_samples": os_, "shift_windows": ow, "r": r}
    return {
        "n_windows_raw": int(n_raw),
        "n_windows_labels": int(n_lab),
        "count_match": bool(n_raw == n_lab),
        "shift_samples": int(best["shift_samples"]),
        "shift_windows": int(best["shift_windows"]),
        "r_min": float(best["r"].min()),
        "r": [float(x) for x in best["r"]],
        "accepted": bool(best["r"].min() >= threshold and n_raw >= n_lab),
    }
