"""Evaluation never leaks the test participant; the decomposition recovers known inputs."""

import numpy as np
import pytest
from conftest import ZHANG, needs_zhang

from edasynth.decompose import bateman, decompose
from edasynth.models import auprc, bootstrap_participants, loso_predict, sensitivity_at_specificity


def test_loso_predictions_ignore_the_test_participants_labels(cfg):
    rng = np.random.default_rng(0)
    X = rng.standard_normal((600, 5))
    g = np.repeat(np.arange(6), 100)
    y = (X[:, 0] + 0.5 * rng.standard_normal(600) > 1).astype(int)
    p1, _ = loso_predict(X, y, g, "logreg", cfg)
    y2 = y.copy()
    y2[g == 2] = 1 - y2[g == 2]  # change labels of participant 2 only
    p2, _ = loso_predict(X, y2, g, "logreg", cfg)
    assert np.allclose(p1[g == 2], p2[g == 2])
    assert not np.allclose(p1[g != 2], p2[g != 2])


def test_metrics_on_known_cases():
    y = np.array([0, 0, 0, 1, 1])
    assert auprc(y, np.array([0.1, 0.2, 0.3, 0.8, 0.9])) == pytest.approx(1.0)
    assert sensitivity_at_specificity(y, np.array([0.1, 0.2, 0.3, 0.8, 0.9]), 0.95) == 1.0


def test_participant_bootstrap_ci_contains_zero_for_identical_scores():
    rng = np.random.default_rng(1)
    y = (rng.random(500) < 0.1).astype(int)
    s = rng.random(500)
    g = np.repeat(np.arange(10), 50)
    lo, hi, d = bootstrap_participants(y, s, s, g, 100, 0)
    assert lo == hi == 0.0


def test_numpy_decomposition_recovers_phasic_component():
    rng = np.random.default_rng(1)
    fs, n = 8.0, 8 * 300
    p = np.zeros(n)
    idx = rng.choice(np.arange(80, n - 80), 12, replace=False)
    p[idx] = rng.uniform(0.5, 2, 12) * fs
    r = np.convolve(p, bateman(np.arange(0, 40, 1 / fs), 3.0))[:n] / fs
    y = r + 1.0 + 0.002 * rng.standard_normal(n)
    z = (y - y.mean()) / y.std()
    out = decompose(z, 1 / fs, 3.0, backend="numpy", max_iter=8000)
    assert np.corrcoef(out["r"], r)[0, 1] > 0.99
    assert np.all(out["p"] >= 0)


def test_numpy_and_reference_backends_reach_the_same_objective():
    from edasynth.decompose import _reference_module

    try:
        _reference_module()
    except ImportError:
        pytest.skip("reference cvxEDA.py not on the path")
    rng = np.random.default_rng(2)
    n = 8 * 120
    y = np.cumsum(rng.standard_normal(n)) * 0.01
    y = (y - y.mean()) / y.std()
    ref = decompose(y, 1 / 8, 3.0, backend="cvxeda")
    ours = decompose(y, 1 / 8, 3.0, backend="numpy", max_iter=30000, tol=1e-12)
    assert ours["obj"] <= ref["obj"] * (1 + 1e-3)


@needs_zhang
@pytest.mark.slow
def test_reproduces_published_feature_baseline(cfg):
    """EDA-only 0.449 vs EDA + |a| 0.475 (logistic regression, LOSO, stored features)."""
    from edasynth.config import merge
    from edasynth.data import load_utd_labels, load_utd_stored_features
    from edasynth.realdata import compare

    lab = load_utd_labels(ZHANG)
    eda, acc = load_utd_stored_features(ZHANG)
    rc = merge(cfg, {"models": cfg["realdata"]["models"], "analysis": {"bootstrap_n": 200}})
    r = compare(eda, np.hstack([eda, acc[:, :4]]), lab.majority.to_numpy(), lab.group.to_numpy(), rc)
    assert r["auprc_0"] == pytest.approx(0.449, abs=0.001)
    assert r["auprc_1"] == pytest.approx(0.475, abs=0.001)
