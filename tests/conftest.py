"""Shared fixtures. Tests that need downloaded data are skipped when it is absent."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from edasynth.config import load_config  # noqa: E402

ZHANG = ROOT / "data" / "external" / "eda-artifact-detection"
NONEEG = ROOT / "data" / "raw" / "noneeg-1.0.0"

needs_zhang = pytest.mark.skipif(not (ZHANG / "UTDallas").exists(), reason="Zhang et al. repository not downloaded")
needs_noneeg = pytest.mark.skipif(
    not (NONEEG / "Subject1_AccTempEDA.hea").exists(), reason="PhysioNet noneeg not downloaded"
)


@pytest.fixture(scope="session")
def cfg():
    return load_config(ROOT / "configs" / "default.yaml")


@pytest.fixture
def toy_bank():
    """A tiny, fully synthetic bank: 4 'participants' with made-up physiology.

    Used only to exercise code paths; its numbers are not estimates of anything.
    """
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
        motion.append({"pid": f"T{k}", "fs": fs, "acc": acc, "stage": stage, "events": events, "q_acc": 0.01})
        phys.append(
            {
                "pid": f"T{k}",
                "tau0": 2.5 + 0.3 * k,
                "rates": {"rest": {"rate_per_s": 0.03}, "load": {"rate_per_s": 0.06}},
                "areas": list(rng.uniform(0.5, 2.0, 50)),
                "tonic": [{"s0": 2.0, "b": 0.001, "A_s": 0.05, "T_t": 60.0, "stage_type": t} for t in ("rest", "load")],
                "sigma": 0.003,
                "q": 0.001,
            }
        )
    lib = {
        "C1": [{"A": 0.5, "tau_r": 0.5, "tau_f": 1.0}],
        "C2": [{"depth": 0.4, "w": 0.125, "D": 2.0}],
        "C3": [{"delta": 0.3, "w": 0.5}],
        "C5": [{"D": 2.0, "kind": "zero"}],
        "class_shares": {"C1": 0.25, "C2": 0.25, "C3": 0.25, "C5": 0.25},
        "G": 0.2,
    }
    return {"motion": motion, "phys": phys, "library": lib, "source": "toy"}
