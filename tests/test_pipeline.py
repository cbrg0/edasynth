"""End-to-end on the toy bank: data set, evaluation, controls, detectors."""

import numpy as np

from edasynth.config import merge
from edasynth.detectors import chen_swt, kleckner, kocielnik, pn_score, slope_score
from edasynth.experiment import evaluate, make_dataset, required_reps
from edasynth.generator import Condition


def _fast(cfg):
    return merge(cfg, {"models": {"logreg": {"C_selection": "fixed"}}})


def test_dataset_shapes(toy_bank, cfg):
    ds = make_dataset(toy_bank, Condition(c=1.0, prevalence=0.05), cfg, seed=0)
    n = ds["y"].size
    assert ds["F0"].shape == (n, 24) and ds["FA"].shape == (n, 96)
    assert ds["FA4"].shape == (n, 4) and ds["FB"].shape == (n, 4)
    assert set(np.unique(ds["groups"])) == {0, 1, 2, 3}


def test_evaluation_rows_have_all_deltas(toy_bank, cfg):
    ds = make_dataset(toy_bank, Condition(c=1.0, prevalence=0.05), cfg, seed=1)
    row = evaluate(ds, _fast(cfg), ("logreg",), 0)[0]
    for k in ("dA", "dA4", "dB", "dB_minus_dA4", "lost_motion_F0", "sens95_F0+FB"):
        assert k in row


def test_acc_shift_breaks_alignment(toy_bank, cfg):
    a = make_dataset(toy_bank, Condition(c=1.0, prevalence=0.05), cfg, seed=2)
    b = make_dataset(toy_bank, Condition(c=1.0, prevalence=0.05), cfg, seed=2, acc_shift=True)
    assert np.array_equal(a["y"], b["y"])
    assert not np.allclose(a["FA4"], b["FA4"])


def test_required_reps_formula():
    assert required_reps(0.01, 0.005, 50) == 16
    assert required_reps(1.0, 0.005, 50) == 50


def test_detectors_run_and_flag_a_dropout():
    fs = 8.0
    n = int(fs * 900)
    y = 3 + 0.01 * np.sin(np.arange(n) / 50)
    y[4000:4020] = 0.0

    cfg = {
        "detectors": {
            "kleckner": {
                "fir_seconds": 64,
                "cutoff_hz": 0.35,
                "level_min": 0.05,
                "level_max": 60,
                "slope_max": 10,
                "transition_s": 5,
            },
            "kocielnik": {"rise_per_s": 0.2, "fall_per_s": 0.1},
            "chen": {"levels": 8, "delta": 0.01, "window_s": 400},
        }
    }
    w = 4000 // 40
    assert kleckner(y, fs, cfg)[w] == 1
    assert kocielnik(y, fs, cfg)[w] == 1
    assert chen_swt(y, fs, cfg)[w] == 1
    assert slope_score(y, fs)[w] == slope_score(y, fs).max()
    assert pn_score(y, fs).shape == (n // 40,)


def test_resumable_results_keep_columns_aligned(tmp_path, toy_bank, cfg):
    import pandas as pd

    from edasynth.cli import _run_jobs
    from edasynth.generator import Condition

    small = merge(cfg, {"models": {"logreg": {"C_selection": "fixed"}}})
    out = tmp_path / "m.csv"
    a = Condition(c=0.0, classes=("C1",))
    b = Condition(c=1.0, classes=("C1", "C2"))
    _run_jobs(small, toy_bank, [(a, 0)], 1, out, ("logreg",))
    _run_jobs(small, toy_bank, [(a, 0), (b, 0)], 1, out, ("logreg",))
    df = pd.read_csv(out)
    assert len(df) == 2
    assert set(df["c"]) == {0.0, 1.0}
    assert set(df["classes"]) == {"C1", "C1+C2"}
