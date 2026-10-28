"""G1: the signal model does exactly what the formulas say."""

import numpy as np
import pytest
from scipy import stats

from edasynth.decompose import bateman, hMax, tPeak
from edasynth.generator import (
    Condition,
    Generator,
    assignDonors,
    blockMean,
    derangement,
    stepShape,
    trapezoid,
    windowLabels,
)


def testBatemanPeakTimeMatchesFormula():
    tau0, tau1 = 3.0, 0.7
    assert tPeak(tau0, tau1) == pytest.approx(1.329, abs=1e-3)
    t = np.arange(0, 20, 1 / 1000)
    assert t[np.argmax(bateman(t, tau0, tau1))] == pytest.approx(tPeak(tau0, tau1), abs=1e-3)
    assert hMax(tau0, tau1) == pytest.approx(bateman(np.array([tPeak(tau0, tau1)]), tau0, tau1)[0])


@pytest.mark.parametrize("tau0,expected", [(2.0, 0.37), (3.0, 0.49), (4.0, 0.57)])
def testBatemanPeakHeight(tau0, expected):
    assert hMax(tau0) == pytest.approx(expected, abs=0.01)


def testResponseCountIsPoisson(toyBank, cfg):
    """Over a long record the number of spontaneous responses lies in the 95 % Poisson interval."""
    gen = Generator(toyBank["library"], cfg)
    m = dict(toyBank["motion"][0])
    n = int(8 * 3600 * 3)
    m["acc"] = np.zeros((n, 3))
    m["stage"] = np.array(["relax1"] * n, dtype=object)
    m["events"] = []
    p = toyBank["phys"][1]
    counts = []
    for seed in range(5):
        cleanSignal, resp = gen.clean(m, p, np.random.default_rng(seed), 0.0)
        counts.append(len(resp))
    lam = p["rates"]["rest"]["ratePerS"] * n / 8
    lo, hi = stats.poisson.ppf([0.0005, 0.9995], lam)
    assert all(lo <= c <= hi for c in counts)


def testTemplatesHaveRequestedHeightAndDuration():
    fs = 64.0
    a = trapezoid(int(fs * 20), fs, 5.0, 0.5, 2.0, 1.0, 0.7)
    assert a.max() == pytest.approx(0.7)
    support = np.where(a > 0)[0] / fs
    assert support[0] == pytest.approx(5.0, abs=2 / fs)
    assert support[-1] == pytest.approx(5.0 + 0.5 + 2.0 + 1.0, abs=2 / fs)
    s = stepShape(int(fs * 20), fs, 3.0, 1.0, -0.4)
    assert s[-1] == pytest.approx(-0.4)
    assert s[int(2.9 * fs)] == 0.0


def testBlockMeanDownsamplesByEight():
    x = np.arange(64.0)
    assert np.allclose(blockMean(x), [3.5 + 8 * k for k in range(8)])


def testWindowLabelRules():
    m = np.zeros(80, bool)
    m[3] = True  # one corrupted sample in window 0
    m[40:65] = True  # 25 of 40 in window 1
    assert list(windowLabels(m, 2, "any")) == [1, 1]
    assert list(windowLabels(m, 2, "majority")) == [0, 1]


def testDerangementHasNoFixedPoint():
    rng = np.random.default_rng(1)
    for iteration in range(200):
        p = derangement(20, rng)
        assert not np.any(p == np.arange(20))
        assert sorted(p) == list(range(20))


def testGeneratedParticipantUsesDifferentDonors(toyBank, cfg):
    from edasynth.experiment import makeDataset

    ds = makeDataset(toyBank, Condition(c=0.5, prevalence=0.05), cfg, seed=3)
    for r in ds["records"]:
        assert r.meta["motionDonor"] != r.meta["physiologyDonor"]


def testPrevalenceAndLabelsFollowTheMask(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    cond = Condition(c=0.0, classes=("C1", "C2", "C5"), prevalence=0.05)
    rec = gen.generate("x", toyBank["motion"][0], toyBank["phys"][1], cond, np.random.default_rng(2))
    assert rec.labels.mean() >= 0.05
    assert np.array_equal(rec.labels, windowLabels(rec.mask, rec.labels.size, "any"))
    diff = np.abs(rec.y - rec.y0)
    q = toyBank["phys"][1]["q"]
    cleanWindows = np.where(rec.labels == 0)[0]
    for w in cleanWindows[:50]:
        assert diff[w * 40 : (w + 1) * 40].max() < cfg["generator"]["labelThreshold"] + 2 * q


def testCZeroMeansNoArtifactIsCausedByMotion(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    rec = gen.generate(
        "x", toyBank["motion"][0], toyBank["phys"][1], Condition(c=0.0, prevalence=0.05), np.random.default_rng(4)
    )
    assert rec.artifacts and not any(a["causedByMotion"] for a in rec.artifacts)
    rec = gen.generate(
        "x",
        toyBank["motion"][0],
        toyBank["phys"][1],
        Condition(c=1.0, deltaMax=1.0, prevalence=0.05),
        np.random.default_rng(4),
    )
    assert all(a["causedByMotion"] for a in rec.artifacts)
    for a in rec.artifacts:
        assert 0 <= a["t0"] - a["eventOnset"] / 8.0 <= 1.0 + 1e-9


def testKappaScalesAmplitudeWithMotionIntensity(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    rec = gen.generate(
        "x",
        toyBank["motion"][0],
        toyBank["phys"][1],
        Condition(c=1.0, kappa=1.0, classes=("C1",), prevalence=0.05),
        np.random.default_rng(5),
    )
    ev = {e["onset"]: e["intensity"] for e in toyBank["motion"][0]["events"]}
    med = np.median(list(ev.values()))
    for a in rec.artifacts:
        assert a["scale"] == pytest.approx(ev[a["eventOnset"]] / med)


def testLevelShiftLabelsOnlyTheTransition(toyBank, cfg):
    """C3 leaves a lasting offset, but only the transition windows are labelled."""
    gen = Generator(toyBank["library"], cfg)
    rec = gen.generate(
        "x",
        toyBank["motion"][0],
        toyBank["phys"][1],
        Condition(c=0.0, classes=("C3",), prevalence=0.01),
        np.random.default_rng(6),
    )
    offset = np.abs(rec.y - rec.y0) >= cfg["generator"]["labelThreshold"]
    assert offset.mean() > rec.mask.mean()  # the offset persists after the labelled part
    assert rec.labels.mean() < 0.05


def testNoiseIsMadeOfWholeResidualBlocks(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    p = dict(toyBank["phys"][0])
    p["residual"] = np.arange(3 * 40, dtype=float).reshape(3, 40)
    e = gen.noise(p, 130, np.random.default_rng(0))
    assert e.size == 130
    for k in range(3):  # every full block is one of the stored blocks, in order
        assert any(np.array_equal(e[40 * k : 40 * k + 40], b) for b in p["residual"])
    p["residual"] = np.zeros((0, 40))
    with pytest.raises(ValueError):
        gen.noise(p, 10, np.random.default_rng(0))


def testArtifactSizeScalesWithTheLocalLevel(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    n, fs = 8 * 60, 8.0
    acc = np.zeros((n, 3))
    heights = []
    for level in (1.0, 3.0, 0.01):
        y0 = np.full(n, level)
        c1 = {"relA": 0.2, "tauR": 0.5, "tauF": 1.0}
        a, sup = gen.artifactSignal("C1", c1, 30.0, 1.0, n, fs, y0, acc, None, "V1")
        heights.append(a.max())
        d, sup = gen.artifactSignal("C3", {"relDelta": -0.1, "w": 0.5}, 30.0, 1.0, n, fs, y0, acc, None, "V1")
        assert d[-1] == pytest.approx(-0.1 * max(level, 0.05))
    assert heights[0] == pytest.approx(0.2, rel=0.1)  # 8-Hz averaging lowers a sharp peak slightly
    assert heights[1] == pytest.approx(3 * heights[0])
    assert heights[2] == pytest.approx(0.05 * heights[0])  # level below the floor


def testInvalidPhysiologyDonorsAreNeverUsed(toyBank):
    phys = [dict(p) for p in toyBank["phys"]]
    phys[2]["validDonor"] = False
    rng = np.random.default_rng(3)
    for iteration in range(200):
        perm = assignDonors(toyBank["motion"], phys, rng)
        assert 2 not in perm
        assert not np.any(perm == np.arange(4))
        counts = np.bincount(perm, minlength=4)[[0, 1, 3]]
        assert counts.max() - counts.min() <= 1
    perm = assignDonors(toyBank["motion"], toyBank["phys"], np.random.default_rng(3))
    assert sorted(perm) == [0, 1, 2, 3]  # all valid: a derangement


def testC4ControlGainGivesTheTargetSize(toyBank, cfg):
    from edasynth.generator import c4ControlGain

    fs, n = 8.0, 8 * 120
    t = np.arange(n) / fs
    acc = np.column_stack([np.zeros(n), np.zeros(n), np.ones(n)])
    events = []
    for s in range(80, n - 80, 160):  # identical bursts of a 2-Hz oscillation
        acc[s : s + 40, 2] += 0.3 * np.sin(2 * np.pi * 2.0 * t[s : s + 40])
        events.append({"onset": s, "end": s + 40, "intensity": 1.0})
    motion = [{"pid": "M", "fs": fs, "acc": acc, "events": events}]
    G = c4ControlGain(motion, cfg["estimation"]["c4"]["smoothS"], 0.1)
    lib = dict(toyBank["library"], Gcontrol=G)
    gen = Generator(lib, cfg)
    out, sup = gen.artifactSignal("C4", {}, 0.0, 1.0, n, fs, np.ones(n), acc, events[1], "V1")
    e = events[1]
    assert np.median(np.abs(out[e["onset"] : e["end"]])) == pytest.approx(0.1, rel=1e-6)


def testTransplantedEpisodeKeepsItsMeasuredCourse(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    n, fs = 8 * 60, 8.0
    y0 = np.full(n, 2.0)
    acc = np.zeros((n, 3))
    trace = [0.0, 0.1, 0.3, 0.2, 0.05, 0.0]
    a, sup = gen.artifactSignal(
        "C1", {"relA": 0.3, "tauR": 0.1, "tauF": 0.1, "trace": trace}, 10.0, 1.0, n, fs, y0, acc, None, "V1"
    )
    assert np.allclose(a[80:86], 2.0 * np.array(trace)) and np.all(a[86:] == 0)
    assert list(np.where(sup)[0]) == [81, 82, 83, 84]  # |piece| >= 0.05 uS inside the episode
    step = [0.0, -0.1, -0.2, -0.2]
    d, sup = gen.artifactSignal(
        "C3", {"relDelta": -0.2, "w": 0.2, "trace": step}, 10.0, 1.0, n, fs, y0, acc, None, "V1"
    )
    assert np.allclose(d[84:], -0.4) and not sup[84:].any()  # lasting offset, not labelled
    deep, sup = gen.artifactSignal(
        "C2", {"relDepth": 0.9, "w": 0.1, "D": 1.0, "trace": [-0.9] * 8}, 10.0, 5.0, n, fs, y0, acc, None, "V1"
    )
    assert deep.min() == pytest.approx(-cfg["generator"]["rhoMax"] * 2.0)  # bounded after movement scaling


def testTonicFollowsTheDonorsOwnLevels(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    fs = 8.0
    n = int(fs * 600)
    m = dict(toyBank["motion"][0], acc=np.zeros((n, 3)), events=[])
    m["stage"] = np.array(["relax1"] * (n // 2) + ["physical"] * (n - n // 2), dtype=object)
    p = dict(toyBank["phys"][0], rates={"rest": {"ratePerS": 0.0}, "load": {"ratePerS": 0.0}})
    p["residual"] = np.zeros((4, 40))
    p["tonic"] = [  # end levels 1.0 + 0.002 * 400 = 1.8 and 3.0 - 0.001 * 1000 = 2.0
        {"s0": 1.0, "b": 0.002, "AS": 0.0, "TT": 60.0, "stageType": "rest", "stage": "relax1", "durationS": 400.0},
        {"s0": 3.0, "b": -0.001, "AS": 0.0, "TT": 60.0, "stageType": "load", "stage": "physical", "durationS": 1000.0},
    ]
    y0, resp = gen.clean(m, p, np.random.default_rng(0), 0.0)
    half = n // 2
    assert y0[0] == pytest.approx(1.0, abs=0.01)
    assert y0[half - 1] == pytest.approx(1.8, abs=0.01)  # end of the first stage, whatever its length here
    assert y0[-1] == pytest.approx(2.0, abs=0.01)  # no jump between stages, ends at the donor's level
    assert np.max(np.abs(np.diff(y0))) < 0.01


def testLevelShiftsReplaceEachOtherInsteadOfAddingUp(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    n = 8 * 120
    y0 = np.full(n, 2.0)
    first = {"i0": 80, "course": np.array([0.0, 0.1, 0.2]), "new": 0.2, "lvl": 2.0}
    second = {"i0": 400, "course": np.array([0.0, -0.05, -0.1]), "new": -0.1, "lvl": 2.0}
    off, sup = gen.contactCourse([second, first], n, y0, "V1")  # order of insertion does not matter
    assert off[200] == pytest.approx(0.4)  # 2 uS x 0.2 after the first shift
    assert off[-1] == pytest.approx(-0.2)  # the second shift replaces the state: 2 x -0.1, not 0.4 - 0.2
    assert not sup[83:400].any() and not sup[403:].any()  # lasting offsets are not labelled
    assert sup[81:83].all() and sup[400:403].any()


def testManyLevelShiftsKeepTheSignalInRange(toyBank, cfg):
    lib = dict(toyBank["library"], C3=[{"relDelta": 0.5, "w": 0.2, "trace": [0.0, 0.25, 0.5]}])
    gen = Generator(lib, cfg)
    rec = gen.generate(
        "x",
        toyBank["motion"][0],
        toyBank["phys"][1],
        Condition(c=0.0, classes=("C3",), prevalence=0.3),
        np.random.default_rng(8),
    )
    # 30 % of windows need many shifts. Added up they would reach several times the
    # level; with replacement the offset stays near +50 % of the level at each shift
    # (a little more while one shift blends into the next, and where y0 has fallen since)
    assert len(rec.artifacts) > 10
    assert np.max(rec.y / rec.y0) < 2.0
    assert np.max(rec.y) <= cfg["generator"]["saturationLevel"]


def testVariedResponseShapesKeepTheHeights(toyBank, cfg):
    gen = Generator(toyBank["library"], cfg)
    m = dict(toyBank["motion"][0])
    p = dict(toyBank["phys"][1], residual=np.zeros((4, 40)))
    a, ra = gen.clean(m, p, np.random.default_rng(3), 0.0)
    gen.responseJitter = cfg["validation"]["g2Jitter"]
    b, rb = gen.clean(m, p, np.random.default_rng(3), 0.0)
    assert len(ra) == len(rb) and not np.allclose(a, b)  # same responses, different shapes
    assert abs(a.max() - b.max()) < 0.2 * (a.max() - a.min())
