"""Estimation (docs/METHOD.md section 3) runs on a labelled recording and recovers what was planted."""

import numpy as np
import pytest

from edasynth.config import merge
from edasynth.data import Recording
from edasynth.estimation import (
    artifactEpisodes,
    cleanEvents,
    cleanWindows,
    decomposeParticipant,
    donorLevelCheck,
    estimateParticipant,
    motionEvents,
    noiseSigma,
    quantisationStep,
    residualBlocks,
    templateParameters,
)
from edasynth.generator import Condition, Generator


def fastCfg(cfg):
    return merge(cfg, {"estimation": {"cvxeda": {"backend": "numpy", "tau0Min": 2.5, "tau0Max": 3.5, "tau0Step": 0.5}}})


def testQuantisationStep():
    x = np.round(np.cumsum(np.random.default_rng(0).integers(-2, 3, 1000)) * 0.003, 6)
    assert quantisationStep(x) == pytest.approx(0.003)


def testCleanWindowsNeedCleanNeighbours():
    lab = np.zeros((6, 3), int)
    lab[2, 0] = 1
    r = Recording("x", 8.0, np.zeros(240), np.zeros((240, 3)), labels3=lab)
    assert list(cleanWindows(r)) == [True, False, False, False, True, True]


def testEstimationEndToEndOnASyntheticRecording(toyBank, cfg):
    cfg = fastCfg(cfg)
    gen = Generator(toyBank["library"], cfg)
    m, p = toyBank["motion"][0], toyBank["phys"][1]  # tau0 = 2.8 s
    rec = gen.generate("x", m, p, Condition(c=1.0, classes=("C1", "C2"), prevalence=0.05), np.random.default_rng(1))
    lab3 = np.repeat(rec.labels[:, None], 3, axis=1)
    r = Recording("x", rec.fs, rec.y, rec.acc, stage=rec.stage, labels3=lab3)
    est = estimateParticipant(r, cfg)
    assert est["tau0"] in (2.5, 3.0, 3.5)
    assert abs(est["tau0"] - p["tau0"]) <= 0.5
    assert est["sigma"] == pytest.approx(p["sigma"], rel=0.5)
    assert est["residual"].shape[1] == 40 and est["residualSource"] == "cleanRest"
    assert est["validDonor"]
    assert len(est["episodes"]) > 0
    assert {e["cls"] for e in est["episodes"]} <= {"C1", "C2", "C3", "C4", "C5", "unclassified"}


def testMotionEventsFoundWherePlanted(cfg):
    rng = np.random.default_rng(2)
    n = 8 * 600
    acc = np.tile([0.0, 0.0, 1.0], (n, 1)) + 0.002 * rng.standard_normal((n, 3))
    stage = np.array(["relax1"] * (n // 2) + ["physical"] * (n - n // 2), dtype=object)
    planted = list(range(n // 2 + 200, n - 100, 300))
    for s in planted:
        acc[s : s + 24] += 0.4 * rng.standard_normal((24, 3))
    r = Recording("x", 8.0, np.ones(n), acc, stage=stage)
    ev, unused = motionEvents(r, cfg)
    found = [e["onset"] for e in ev]
    assert all(any(abs(f - p) <= 4 for f in found) for p in planted)
    assert len(ev) == len(planted)


def testNoiseSigmaRecoversWhiteNoise():
    rng = np.random.default_rng(3)
    n = 40 * 200
    y = 2 + 0.01 * rng.standard_normal(n)
    r = Recording(
        "x", 8.0, y, np.zeros((n, 3)), stage=np.array(["relax1"] * n, dtype=object), labels3=np.zeros((200, 3), int)
    )
    assert noiseSigma(r, np.ones(200, bool)) == pytest.approx(0.01, rel=0.05)


def testEpisodeClassificationRules(cfg):
    n = 40 * 60
    y = np.full(n, 3.0)
    y[400:440] -= 0.5  # contact drop that returns
    y[1600:] += 0.4  # lasting level shift
    lab = np.zeros((60, 3), int)
    lab[10, :] = 1
    lab[40, :] = 1
    r = Recording("x", 8.0, y, np.zeros((n, 3)), stage=np.array(["relax1"] * n, dtype=object), labels3=lab)
    eps = artifactEpisodes(r, 0.001, cfg)
    assert [e["cls"] for e in eps] == ["C2", "C3"]


def testNoiseSigmaFallsBackWhenNoRestWindows():
    rng = np.random.default_rng(5)
    n = 8 * 300
    eda = 1.0 + 0.01 * rng.standard_normal(n)
    r = Recording("w", 8.0, eda, np.zeros((n, 3)), stage=np.full(n, "walk", dtype=object))
    clean = np.ones(n // 40, bool)
    sigma, src = noiseSigma(r, clean, withSource=True)
    assert src == "cleanAllFallback"
    assert 0.007 < sigma < 0.013


@pytest.mark.parametrize("tauTrue", [2.5, 4.0])
def testTau0SelectionRecoversKnownValue(cfg, tauTrue):
    from edasynth.decompose import bateman

    fs, n = 8.0, 8 * 600
    t = np.arange(n) / fs
    rng = np.random.default_rng(1)
    y = 2 + 0.001 * t + 0.003 * rng.standard_normal(n)
    for o in rng.uniform(0, 590, 15):
        y += rng.uniform(0.05, 0.5) * bateman(t - o, tauTrue)
    y = 0.001 * np.round(y / 0.001)
    c = merge(cfg, {"estimation": {"cvxeda": {"backend": "numpy", "tau0Min": 2.0, "tau0Max": 4.0, "tau0Step": 0.5}}})
    r = Recording("x", fs, y, np.zeros((n, 3)), stage=np.full(n, "relax1", dtype=object))
    dec = decomposeParticipant(r, np.ones(n // 40, bool), c)
    assert abs(dec["tau0"] - tauTrue) <= 0.5


def testResidualBlocksComeFromCleanRestWindowsAndAreCentred():
    n = 40 * 6
    stage = np.array(["relax1"] * 120 + ["physical"] * 120, dtype=object)
    r = Recording("x", 8.0, np.ones(n), np.zeros((n, 3)), stage=stage, labels3=np.zeros((6, 3), int))
    resid = np.repeat(np.arange(6.0), 40) + np.tile([0.01, -0.01], 120)
    clean = np.array([True, False, True, True, True, True])
    out = residualBlocks(r, resid, clean, 5.0)
    assert out["source"] == "cleanRest"
    assert out["blocks"].shape == (2, 40)  # windows 0 and 2; window 1 is not clean, 3-5 are load
    assert np.allclose(out["blocks"].mean(1), 0.0)  # window offsets 0 and 2 removed
    assert np.allclose(np.abs(out["blocks"]), 0.01)
    out = residualBlocks(r, resid, np.array([False, False, False, True, True, False]), 5.0)
    assert out["source"] == "cleanAllFallback" and out["blocks"].shape == (2, 40)
    with pytest.raises(ValueError):
        residualBlocks(r, resid, clean, 3.0)  # 24 samples do not divide a 40-sample window


def testTemplateSizesAreRelativeToTheLevelBeforeTheEpisode():
    ep = {"riseS": 0.4, "returnS": 0.8, "durationS": 5.0, "amp": -0.6, "shift": 0.3, "level": 3.0}
    assert templateParameters({**ep, "cls": "C1"}, 0.05)["relA"] == pytest.approx(0.2)
    assert templateParameters({**ep, "cls": "C2"}, 0.05)["relDepth"] == pytest.approx(0.2)
    assert templateParameters({**ep, "cls": "C3"}, 0.05)["relDelta"] == pytest.approx(0.1)
    low = templateParameters({**ep, "cls": "C3", "level": 0.01}, 0.05)
    assert low["relDelta"] == pytest.approx(0.3 / 0.05)  # the floor applies below 0.05 uS


def testEpisodesStoreTheLevelBeforeAndAfter(cfg):
    fs, n = 8.0, 8 * 60
    y = np.full(n, 2.0)
    y[200:240] = 2.5
    y[240:] = 2.2
    lab = np.zeros((n // 40, 3), int)
    lab[5] = 1
    eps = artifactEpisodes(Recording("x", fs, y, np.zeros((n, 3)), labels3=lab), 0.001, cfg)
    assert len(eps) == 1
    assert eps[0]["level"] == pytest.approx(2.0) and eps[0]["levelAfter"] == pytest.approx(2.2)


def testDonorLevelCheckFollowsTheMinimumLevelRule(cfg):
    ok = donorLevelCheck(np.r_[np.full(60, 0.04), np.full(40, 1.0)], cfg)
    assert ok["shareBelowMinLevel"] == pytest.approx(0.6) and not ok["validDonor"]
    ok = donorLevelCheck(np.r_[np.full(50, 0.04), np.full(50, 1.0)], cfg)
    assert ok["validDonor"]  # exactly half below is still accepted


def testResponseSizesComeOnlyFromCleanWindows():
    clean = np.array([True, False, True, True])
    ev = [
        {"onset": 5, "peak": 20},  # all in window 0: kept
        {"onset": 30, "peak": 45},  # peak in the artifact window 1: dropped
        {"onset": 50, "peak": 60},  # onset in window 1: dropped
        {"onset": 90, "peak": 130},  # windows 2 and 3: kept
    ]
    assert [e["onset"] for e in cleanEvents(ev, clean)] == [5, 90]
