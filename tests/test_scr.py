import numpy as np

from edasynth.decompose import bateman, h_max
from edasynth.estimation import scr_events


def test_scr_events_recovers_isolated_responses():
    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 120)) / fs
    onsets, amps = [10.0, 40.0, 75.0], [0.3, 0.1, 0.6]
    r = np.zeros_like(t)
    for o, A in zip(onsets, amps, strict=True):
        r += A / h_max(tau0) * bateman(t - o, tau0)
    ev = scr_events(r, fs, tau0, 0.7, 0.05)
    assert len(ev) == 3
    for e, o, A in zip(ev, onsets, amps, strict=True):
        assert abs(e["onset"] / fs - o) <= 0.25
        assert abs(e["amp"] - A) / A < 0.02
        assert abs(e["area"] * h_max(tau0) - e["amp"]) < 1e-12


def test_scr_events_ignores_small_bumps():
    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 60)) / fs
    r = 0.02 / h_max(tau0) * bateman(t - 10, tau0)
    assert scr_events(r, fs, tau0, 0.7, 0.05) == []


def test_events_to_phasic_rebuilds_the_responses():
    from edasynth.estimation import events_to_phasic

    fs, tau0 = 8.0, 3.0
    t = np.arange(int(fs * 120)) / fs
    r = sum(A / h_max(tau0) * bateman(t - o, tau0) for o, A in [(10.0, 0.3), (60.0, 0.2)])
    ev = scr_events(r, fs, tau0, 0.7, 0.05)
    rebuilt = events_to_phasic(ev, t.size, fs, tau0, 0.7)
    assert np.max(np.abs(rebuilt - r)) < 0.01
