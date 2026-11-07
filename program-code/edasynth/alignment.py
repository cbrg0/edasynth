"""This file belongs to the first step, preparing the data. It checks that every expert
label lies on the right piece of the recording before any analysis uses the labels.

The check is needed because Zhang et al. (2017) published their labels separately from
the PhysioNet recordings, and neither source states which label window matches which
samples. Together with the labels, Zhang et al. published a table of simple features of
each window. The file recomputes four of these features from the raw recording with the
formulas of their script and correlates them with the published values. If the labels
were shifted, the correlation would drop. A participant is accepted when all four
correlations reach 0.999; this threshold is a decision of this study and was set so high
because a correct alignment should reproduce the published numbers almost exactly.

The outcome is the list of accepted participants with the shift between labels and
signal. All 20 participants were accepted without a shift.
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
