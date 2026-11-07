"""This file prepares what the tests share, such as a tiny invented set of participants.
Tests that need the downloaded data are skipped when the data are missing.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "program-code"))

from edasynth.config import loadConfig  # noqa: E402

ZHANG = ROOT / "downloaded-data" / "external" / "eda-artifact-detection"
NONEEG = ROOT / "downloaded-data" / "raw" / "noneeg-1.0.0"

needsZhang = pytest.mark.skipif(not (ZHANG / "UTDallas").exists(), reason="Zhang et al. repository not downloaded")
needsNoneeg = pytest.mark.skipif(
    not (NONEEG / "Subject1_AccTempEDA.hea").exists(), reason="PhysioNet noneeg not downloaded"
)


@pytest.fixture(scope="session")
def cfg():
    return loadConfig(ROOT / "settings" / "default.yaml")


@pytest.fixture
def toyBank():
    rng = np.random.default_rng(0)
    fs = 8.0
    motion, phys = [], []
    for k in range(4):
        n = int(fs * 600)
        acc = np.column_stack([np.zeros(n), np.zeros(n), np.ones(n)]) + 0.01 * rng.standard_normal((n, 3))
        events = []
        for s in range(400, n - 200, 600):
            acc[s : s + 40] += 0.5 * rng.standard_normal((40, 3))
            events.append({"onset": s, "end": s + 40, "intensity": float(5 + rng.random())})
        stage = np.array(["relax1"] * (n // 2) + ["physical"] * (n - n // 2), dtype=object)
        motion.append({"pid": f"T{k}", "fs": fs, "acc": acc, "stage": stage, "events": events, "qAcc": 0.01})
        phys.append(
            {
                "pid": f"T{k}",
                "tau0": 2.5 + 0.3 * k,
                "rates": {"rest": {"ratePerS": 0.03}, "load": {"ratePerS": 0.06}},
                "areas": list(rng.uniform(0.5, 2.0, 50)),
                "tonic": [
                    {"s0": 2.0, "b": 0.001, "AS": 0.05, "TT": 60.0, "stageType": t, "stage": st, "durationS": 300.0}
                    for t, st in (("rest", "relax1"), ("load", "physical"))
                ],
                "sigma": 0.003,
                "residual": 0.003 * rng.standard_normal((60, 40)),
                "validDonor": True,
                "q": 0.001,
            }
        )
    lib = {
        "C1": [{"relA": 0.25, "tauR": 0.5, "tauF": 1.0}],
        "C2": [{"relDepth": 0.2, "w": 0.125, "D": 2.0}],
        "C3": [{"relDelta": 0.15, "w": 0.5}],
        "C5": [{"D": 2.0, "kind": "zero"}],
        "classShares": {"C1": 0.25, "C2": 0.25, "C3": 0.25, "C5": 0.25},
        "G": 0.2,
        "Gcontrol": 0.2,
    }
    return {"motion": motion, "phys": phys, "library": lib, "source": "toy", "bankId": "toy000000000"}
