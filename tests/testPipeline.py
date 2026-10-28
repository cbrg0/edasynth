"""End-to-end on the toy bank: data set, evaluation, controls, detectors."""

import numpy as np
import pytest

from edasynth.config import merge
from edasynth.detectors import chenSwt, kleckner, kocielnik, pnScore, slopeScore
from edasynth.experiment import evaluate, makeDataset, requiredReps
from edasynth.generator import Condition


def fastConfig(cfg):
    return merge(cfg, {"models": {"logreg": {"CSelection": "fixed"}}})


def testDatasetShapes(toyBank, cfg):
    ds = makeDataset(toyBank, Condition(c=1.0, prevalence=0.05), cfg, seed=0)
    n = ds["y"].size
    assert ds["F0"].shape == (n, 24) and ds["FA"].shape == (n, 96)
    assert ds["FA4"].shape == (n, 4) and ds["FB"].shape == (n, 4)
    assert set(np.unique(ds["groups"])) == {0, 1, 2, 3}


def testEvaluationRowsHaveAllDeltas(toyBank, cfg):
    ds = makeDataset(toyBank, Condition(c=1.0, prevalence=0.05), cfg, seed=1)
    row = evaluate(ds, fastConfig(cfg), ("logreg",), 0)[0]
    for k in ("dA", "dA4", "dB", "dBMinusDA4", "lostMotionF0", "sens95F0+FB"):
        assert k in row


def testAccShiftBreaksAlignment(toyBank, cfg):
    a = makeDataset(toyBank, Condition(c=1.0, prevalence=0.05), cfg, seed=2)
    b = makeDataset(toyBank, Condition(c=1.0, prevalence=0.05), cfg, seed=2, accShift=True)
    assert np.array_equal(a["y"], b["y"])
    assert not np.allclose(a["FA4"], b["FA4"])


def testRequiredRepsFormula():
    assert requiredReps(0.01, 0.005, 50) == 16
    assert requiredReps(1.0, 0.005, 50) == 50


def testDetectorsRunAndFlagADropout():
    fs = 8.0
    n = int(fs * 900)
    y = 3 + 0.01 * np.sin(np.arange(n) / 50)
    y[4000:4020] = 0.0

    cfg = {
        "detectors": {
            "kleckner": {
                "firSeconds": 64,
                "cutoffHz": 0.35,
                "levelMin": 0.05,
                "levelMax": 60,
                "slopeMax": 10,
                "transitionS": 5,
            },
            "kocielnik": {"risePerS": 0.2, "fallPerS": 0.1},
            "chen": {"levels": 8, "delta": 0.01, "windowS": 400},
        }
    }
    w = 4000 // 40
    assert kleckner(y, fs, cfg)[w] == 1
    assert kocielnik(y, fs, cfg)[w] == 1
    assert chenSwt(y, fs, cfg)[w] == 1
    assert slopeScore(y, fs)[w] == slopeScore(y, fs).max()
    assert pnScore(y, fs).shape == (n // 40,)


def testResumableResultsKeepColumnsAligned(tmp_path, toyBank, cfg):
    import pandas as pd

    from edasynth.cli import runJobs
    from edasynth.generator import Condition

    small = merge(cfg, {"models": {"logreg": {"CSelection": "fixed"}}})
    out = tmp_path / "m.csv"
    a = Condition(c=0.0, classes=("C1",))
    b = Condition(c=1.0, classes=("C1", "C2"))
    runJobs(small, toyBank, [(a, 0)], 1, out, ("logreg",))
    runJobs(small, toyBank, [(a, 0), (b, 0)], 1, out, ("logreg",))
    df = pd.read_csv(out)
    assert len(df) == 2
    assert set(df["c"]) == {0.0, 1.0}
    assert set(df["classes"]) == {"C1", "C1+C2"}


def testResumeRefusesResultsOfAnotherBank(tmp_path, toyBank, cfg):
    from edasynth.cli import runJobs
    from edasynth.generator import Condition

    small = merge(cfg, {"models": {"logreg": {"CSelection": "fixed"}}})
    out = tmp_path / "m.csv"
    a = Condition(c=0.0, classes=("C1",))
    runJobs(small, toyBank, [(a, 0)], 1, out, ("logreg",))
    other = dict(toyBank, bankId="other0000000")
    with pytest.raises(SystemExit):
        runJobs(small, other, [(a, 0), (a, 1)], 1, out, ("logreg",))


def testPowerConditionsSpreadOverTheGain():
    import pandas as pd

    from edasynth.cli import selectPowerConditions

    comp = pd.DataFrame({"c": np.arange(9.0), "dB": np.arange(9.0)[::-1] / 100})
    chosen = selectPowerConditions(comp, 3)
    assert list(chosen.dB) == [0.02, 0.04, 0.06]
    assert len(selectPowerConditions(comp.head(2), 3)) == 2


def testEquivalentGridCellsAreRunOnceAndFilledIn(cfg):
    import pandas as pd

    from edasynth.experiment import canonicalCondition, expandGrid, gridConditions, gridToRun

    run = gridToRun(cfg)
    assert len(gridConditions(cfg)) == 300 and len(run) == 226
    a = Condition(c=0.0, kappa=1.0, deltaMax=5.0, piR=0.5, classes=("C1",))
    assert canonicalCondition(a, cfg).key() == Condition(c=0.0, kappa=0.0, deltaMax=1.0, piR=0.5, classes=("C1",)).key()
    b = Condition(c=0.5, kappa=1.0, deltaMax=5.0, classes=("C5",))
    assert canonicalCondition(b, cfg).kappa == 0.0 and canonicalCondition(b, cfg).deltaMax == 5.0
    m = Condition(c=0.5, kappa=1.0, classes=("C1", "C2", "C3", "C5"))
    assert canonicalCondition(m, cfg) == m
    rows = [{"key": c.key(), "kappa": c.kappa, "deltaMax": c.deltaMax, "dB": 0.01, "rep": 0} for c in run]
    full = expandGrid(pd.DataFrame(rows), cfg)
    assert sorted(full["key"]) == sorted(c.key() for c in gridConditions(cfg))


def testLabelNoiseGridUsesOnlyTheMixtureCells(cfg):
    from edasynth.cli import gridCells

    cells = gridCells(cfg, False, True)
    assert len(cells) == 50
    assert {c.classes for c in cells} == {tuple(cfg["grid"]["compositions"][-1])}
    assert len(gridCells(cfg, False, False)) == 226


def testPowerDoesNotDependOnTheNumberOfJobs(toyBank, cfg):
    from edasynth.experiment import bootstrapPower

    small = merge(cfg, {"models": {"logreg": {"CSelection": "fixed"}}, "analysis": {"bootstrapN": 50}})
    cond = Condition(c=1.0, classes=("C1",), prevalence=0.05)
    a = bootstrapPower(toyBank, cond, small, 3, 2, nJobs=1)
    b = bootstrapPower(toyBank, cond, small, 3, 2, nJobs=2)
    assert a == b


def testRealismGeneratorNeverSeesTheComparedParticipant(toyBank, cfg):
    from edasynth.validation import withoutParticipants

    def est(pid, cls):
        ep = {
            "pid": pid,
            "cls": cls,
            "level": 2.0,
            "riseS": 0.4,
            "returnS": 0.4,
            "durationS": 5.0,
            "amp": 0.3,
            "shift": 0.2,
            "traceDev": [0.0, 0.1],
            "traceStep": [0.0, 0.1],
        }
        return {"pid": pid, "episodes": [ep], "c4": {"G": 0.0}}

    bank = dict(toyBank, estimates=[est("T0", "C1"), est("T1", "C2"), est("T2", "C3"), est("T3", "C1")])
    sub = withoutParticipants(bank, {"T0", "T3"}, cfg)
    assert [m["pid"] for m in sub["motion"]] == ["T1", "T2"]
    assert [p["pid"] for p in sub["phys"]] == ["T1", "T2"]
    assert sub["library"]["C1"] == [] and len(sub["library"]["C2"]) == 1  # T0's and T3's episodes are gone
