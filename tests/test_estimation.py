"""Estimation (step 3) runs on a labelled recording and recovers what was planted."""

import numpy as np
import pytest

from edasynth.config import merge
from edasynth.data import Recording
from edasynth.estimation import (
    artifact_episodes,
    clean_windows,
    decompose_participant,
    estimate_participant,
    motion_events,
    noise_sigma,
    quantisation_step,
)
from edasynth.generator import Condition, Generator


def _fast_cfg(cfg):
    return merge(
        cfg, {"estimation": {"cvxeda": {"backend": "numpy", "tau0_min": 2.5, "tau0_max": 3.5, "tau0_step": 0.5}}}
    )


def test_quantisation_step():
    x = np.round(np.cumsum(np.random.default_rng(0).integers(-2, 3, 1000)) * 0.003, 6)
    assert quantisation_step(x) == pytest.approx(0.003)


def test_clean_windows_need_clean_neighbours():
    lab = np.zeros((6, 3), int)
    lab[2, 0] = 1
    r = Recording("x", 8.0, np.zeros(240), np.zeros((240, 3)), labels3=lab)
    assert list(clean_windows(r)) == [True, False, False, False, True, True]


def test_estimation_end_to_end_on_a_synthetic_recording(toy_bank, cfg):
    cfg = _fast_cfg(cfg)
    gen = Generator(toy_bank["library"], cfg)
    m, p = toy_bank["motion"][0], toy_bank["phys"][1]  # tau0 = 2.8 s
    rec = gen.generate("x", m, p, Condition(c=1.0, classes=("C1", "C2"), prevalence=0.05), np.random.default_rng(1))
    lab3 = np.repeat(rec.labels[:, None], 3, axis=1)
    r = Recording("x", rec.fs, rec.y, rec.acc, stage=rec.stage, labels3=lab3)
    est = estimate_participant(r, cfg)
    assert est["tau0"] in (2.5, 3.0, 3.5)
    assert abs(est["tau0"] - p["tau0"]) <= 0.5
    assert est["sigma"] == pytest.approx(p["sigma"], rel=0.5)
    assert len(est["episodes"]) > 0
    assert {e["cls"] for e in est["episodes"]} <= {"C1", "C2", "C3", "C4", "C5", "unclassified"}


def test_motion_events_found_where_planted(cfg):
    rng = np.random.default_rng(2)
    n = 8 * 600
    acc = np.tile([0.0, 0.0, 1.0], (n, 1)) + 0.002 * rng.standard_normal((n, 3))
    stage = np.array(["relax1"] * (n // 2) + ["physical"] * (n - n // 2), dtype=object)
    planted = list(range(n // 2 + 200, n - 100, 300))
    for s in planted:
        acc[s : s + 24] += 0.4 * rng.standard_normal((24, 3))
    r = Recording("x", 8.0, np.ones(n), acc, stage=stage)
    ev, _ = motion_events(r, cfg)
    found = [e["onset"] for e in ev]
    assert all(any(abs(f - p) <= 4 for f in found) for p in planted)
    assert len(ev) == len(planted)


def test_noise_sigma_recovers_white_noise():
    rng = np.random.default_rng(3)
    n = 40 * 200
    y = 2 + 0.01 * rng.standard_normal(n)
    r = Recording(
        "x", 8.0, y, np.zeros((n, 3)), stage=np.array(["relax1"] * n, dtype=object), labels3=np.zeros((200, 3), int)
    )
    assert noise_sigma(r, np.ones(200, bool)) == pytest.approx(0.01, rel=0.05)


def test_episode_classification_rules(cfg):
    n = 40 * 60
    y = np.full(n, 3.0)
    y[400:440] -= 0.5  # contact drop that returns
    y[1600:] += 0.4  # lasting level shift
    lab = np.zeros((60, 3), int)
    lab[10, :] = 1
    lab[40, :] = 1
    r = Recording("x", 8.0, y, np.zeros((n, 3)), stage=np.array(["relax1"] * n, dtype=object), labels3=lab)
    eps = artifact_episodes(r, 0.001, cfg)
    assert [e["cls"] for e in eps] == ["C2", "C3"]


def test_noise_sigma_falls_back_when_no_rest_windows():
    rng = np.random.default_rng(5)
    n = 8 * 300
    eda = 1.0 + 0.01 * rng.standard_normal(n)
    r = Recording("w", 8.0, eda, np.zeros((n, 3)), stage=np.full(n, "walk", dtype=object))
    clean = np.ones(n // 40, bool)
    sigma, src = noise_sigma(r, clean, with_source=True)
    assert src == "clean_all_fallback"
    assert 0.007 < sigma < 0.013


@pytest.mark.parametrize("tau_true", [2.5, 4.0])
def test_tau0_selection_recovers_known_value(cfg, tau_true):
    from edasynth.decompose import bateman

    fs, n = 8.0, 8 * 600
    t = np.arange(n) / fs
    rng = np.random.default_rng(1)
    y = 2 + 0.001 * t + 0.003 * rng.standard_normal(n)
    for o in rng.uniform(0, 590, 15):
        y += rng.uniform(0.05, 0.5) * bateman(t - o, tau_true)
    y = 0.001 * np.round(y / 0.001)
    c = merge(cfg, {"estimation": {"cvxeda": {"backend": "numpy", "tau0_min": 2.0, "tau0_max": 4.0, "tau0_step": 0.5}}})
    r = Recording("x", fs, y, np.zeros((n, 3)), stage=np.full(n, "relax1", dtype=object))
    dec = decompose_participant(r, np.ones(n // 40, bool), c)
    assert abs(dec["tau0"] - tau_true) <= 0.5
