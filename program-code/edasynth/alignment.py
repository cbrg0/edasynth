"""Linking the expert labels (Zhang et al. repository) to the raw UTD signals.

Labels and signals were published separately and no document states which
window belongs to which piece of signal. Two checks are applied per participant:
1. the number of complete 5-s windows equals the number of labelled windows of
   group g = subject - 1;
2. the four amplitude features recomputed from the raw signal and z-scored within
   the participant (exactly as in Zhang's script) correlate >= 0.999 with the
   first four columns of the published matrix. If not, shifts of 1-10 windows and
   0-39 samples are searched; a participant that still fails is excluded.
Accepted participants are used with the offset found here (``offsetSamples``):
every later step drops that many samples from the start of the raw record.
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
    """Return window counts, best shift and per-feature correlations."""
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
