import numpy as np

from edasynth.decompose import bateman, hMax
from edasynth.estimation import scrEvents


def testScrEventsRecoversIsolatedResponses():
    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 120)) / fs
    onsets, amps = [10.0, 40.0, 75.0], [0.3, 0.1, 0.6]
    r = np.zeros_like(t)
    for o, A in zip(onsets, amps, strict=True):
        r += A / hMax(tau0) * bateman(t - o, tau0)
    ev = scrEvents(r, fs, tau0, 0.7, 0.05)
    assert len(ev) == 3
    for e, o, A in zip(ev, onsets, amps, strict=True):
        assert abs(e["onset"] / fs - o) <= 0.25
        assert abs(e["amp"] - A) / A < 0.02
        assert abs(e["area"] * hMax(tau0) - e["amp"]) < 1e-12


def testScrEventsIgnoresSmallBumps():
    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 60)) / fs
    r = 0.02 / hMax(tau0) * bateman(t - 10, tau0)
    assert scrEvents(r, fs, tau0, 0.7, 0.05) == []


def testEventsToPhasicRebuildsTheResponses():
    from edasynth.estimation import eventsToPhasic

    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 120)) / fs
    r = sum(A / hMax(tau0) * bateman(t - o, tau0) for o, A in [(10.0, 0.3), (60.0, 0.2)])
    ev = scrEvents(r, fs, tau0, 0.7, 0.05)
    rebuilt = eventsToPhasic(ev, t.size, fs, tau0, 0.7)
    assert np.max(np.abs(rebuilt - r)) < 0.01
