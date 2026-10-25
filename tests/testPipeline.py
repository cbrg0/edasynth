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
