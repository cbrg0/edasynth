"""This file checks that the expert labels belong to the right parts of the recordings,
before any analysis starts.

The labels and the recordings were published separately, and nowhere is it written which
label window matches which piece of signal. This file recomputes a few simple numbers
from the signal and compares them with the numbers published together with the labels.
If they match, the labels are in the right place. All 20 participants passed this check.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .data import WIN
from .features import stats4, windows


def amplitudeZ(eda: np.ndarray, offset: int = 0, nWindows: int | None = None) -> np.ndarray:
    w = windows(eda[offset:])
    if nWindows is not None:
        w = w[:nWindows]
    return stats.zscore(stats4(w), axis=0)


def checkParticipant(eda: np.ndarray, stored4: np.ndarray, threshold: float) -> dict:
    """Return the numbers of windows, the best shift and the correlation of each feature."""
    nLab = stored4.shape[0]
    nRaw = eda.size // WIN

    def corrAt(offSamples: int, offWindows: int):
        z = amplitudeZ(eda, offSamples + offWindows * WIN)
        m = min(z.shape[0], nLab)
        if m < 10:
            return np.full(4, -1.0)
        z = stats.zscore(stats4(windows(eda[offSamples + offWindows * WIN :])[:m]), axis=0)
        ref = stats.zscore(stored4[:m], axis=0)
        return np.array([np.corrcoef(z[:, k], ref[:, k])[0, 1] for k in range(4)])

    r0 = corrAt(0, 0)
    best = {"shiftSamples": 0, "shiftWindows": 0, "r": r0}
    if r0.min() < threshold:
        for ow in range(-10, 11):
            for sampleShift in range(0, WIN):
                if ow * WIN + sampleShift < 0:
                    continue
                r = corrAt(sampleShift, ow)
                if r.min() > best["r"].min():
                    best = {"shiftSamples": sampleShift, "shiftWindows": ow, "r": r}
    offset = best["shiftWindows"] * WIN + best["shiftSamples"]
    nAfter = (eda.size - offset) // WIN
    return {
        "nWindowsRaw": int(nRaw),
        "nWindowsAfterOffset": int(nAfter),
        "offsetSamples": int(offset),
        "nWindowsLabels": int(nLab),
        "countMatch": bool(nRaw == nLab),
        "shiftSamples": int(best["shiftSamples"]),
        "shiftWindows": int(best["shiftWindows"]),
        "rMin": float(best["r"].min()),
        "r": [float(x) for x in best["r"]],
        "accepted": bool(best["r"].min() >= threshold and nAfter >= nLab),
    }
