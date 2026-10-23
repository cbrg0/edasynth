"""G1: the signal model does exactly what the formulas say."""

import numpy as np
import pytest
from scipy import stats

from edasynth.decompose import bateman, h_max, t_peak
from edasynth.generator import Condition, Generator, _block_mean, _step, _trapezoid, _window_labels, derangement


def test_bateman_peak_time_matches_formula():
    tau0, tau1 = 3.0, 0.7
    assert t_peak(tau0, tau1) == pytest.approx(1.329, abs=1e-3)
    t = np.arange(0, 20, 1 / 1000)
    assert t[np.argmax(bateman(t, tau0, tau1))] == pytest.approx(t_peak(tau0, tau1), abs=1e-3)
    assert h_max(tau0, tau1) == pytest.approx(bateman(np.array([t_peak(tau0, tau1)]), tau0, tau1)[0])


@pytest.mark.parametrize("tau0,expected", [(2.0, 0.37), (3.0, 0.49), (4.0, 0.57)])
def test_bateman_peak_height(tau0, expected):
    assert h_max(tau0) == pytest.approx(expected, abs=0.01)


def test_response_count_is_poisson(toy_bank, cfg):
    """Over a long record the number of spontaneous responses lies in the 95 % Poisson interval."""
    gen = Generator(toy_bank["library"], cfg)
    m = dict(toy_bank["motion"][0])
    n = int(8 * 3600 * 3)
    m["acc"] = np.zeros((n, 3))
    m["stage"] = np.array(["relax1"] * n, dtype=object)
    m["events"] = []
    p = toy_bank["phys"][1]
    counts = []
    for seed in range(5):
        _, resp = gen.clean(m, p, np.random.default_rng(seed), 0.0)
        counts.append(len(resp))
    lam = p["rates"]["rest"]["rate_per_s"] * n / 8
    lo, hi = stats.poisson.ppf([0.0005, 0.9995], lam)
    assert all(lo <= c <= hi for c in counts)


def test_templates_have_requested_height_and_duration():
    fs = 64.0
    a = _trapezoid(int(fs * 20), fs, 5.0, 0.5, 2.0, 1.0, 0.7)
    assert a.max() == pytest.approx(0.7)
    support = np.where(a > 0)[0] / fs
    assert support[0] == pytest.approx(5.0, abs=2 / fs)
    assert support[-1] == pytest.approx(5.0 + 0.5 + 2.0 + 1.0, abs=2 / fs)
    s = _step(int(fs * 20), fs, 3.0, 1.0, -0.4)
    assert s[-1] == pytest.approx(-0.4)
    assert s[int(2.9 * fs)] == 0.0


def test_block_mean_downsamples_by_eight():
    x = np.arange(64.0)
    assert np.allclose(_block_mean(x), [3.5 + 8 * k for k in range(8)])


def test_window_label_rules():
    m = np.zeros(80, bool)
    m[3] = True  # one corrupted sample in window 0
    m[40:65] = True  # 25 of 40 in window 1
    assert list(_window_labels(m, 2, "any")) == [1, 1]
    assert list(_window_labels(m, 2, "majority")) == [0, 1]


def test_derangement_has_no_fixed_point():
    rng = np.random.default_rng(1)
    for _ in range(200):
        p = derangement(20, rng)
        assert not np.any(p == np.arange(20))
        assert sorted(p) == list(range(20))


def test_generated_participant_uses_different_donors(toy_bank, cfg):
    from edasynth.experiment import make_dataset

    ds = make_dataset(toy_bank, Condition(c=0.5, prevalence=0.05), cfg, seed=3)
    for r in ds["records"]:
        assert r.meta["motion_donor"] != r.meta["physiology_donor"]


def test_prevalence_and_labels_follow_the_mask(toy_bank, cfg):
    gen = Generator(toy_bank["library"], cfg)
    cond = Condition(c=0.0, classes=("C1", "C2", "C5"), prevalence=0.05)
    rec = gen.generate("x", toy_bank["motion"][0], toy_bank["phys"][1], cond, np.random.default_rng(2))
    assert rec.labels.mean() >= 0.05
    assert np.array_equal(rec.labels, _window_labels(rec.mask, rec.labels.size, "any"))
    diff = np.abs(rec.y - rec.y0)
    q = toy_bank["phys"][1]["q"]
    clean_windows = np.where(rec.labels == 0)[0]
    for w in clean_windows[:50]:
        assert diff[w * 40 : (w + 1) * 40].max() < cfg["generator"]["label_threshold"] + 2 * q


def test_c_zero_means_no_artifact_is_caused_by_motion(toy_bank, cfg):
    gen = Generator(toy_bank["library"], cfg)
    rec = gen.generate(
        "x", toy_bank["motion"][0], toy_bank["phys"][1], Condition(c=0.0, prevalence=0.05), np.random.default_rng(4)
    )
    assert rec.artifacts and not any(a["caused_by_motion"] for a in rec.artifacts)
    rec = gen.generate(
        "x",
        toy_bank["motion"][0],
        toy_bank["phys"][1],
        Condition(c=1.0, delta_max=1.0, prevalence=0.05),
        np.random.default_rng(4),
    )
    assert all(a["caused_by_motion"] for a in rec.artifacts)
    for a in rec.artifacts:
        assert 0 <= a["t0"] - a["event_onset"] / 8.0 <= 1.0 + 1e-9


def test_kappa_scales_amplitude_with_motion_intensity(toy_bank, cfg):
    gen = Generator(toy_bank["library"], cfg)
    rec = gen.generate(
        "x",
        toy_bank["motion"][0],
        toy_bank["phys"][1],
        Condition(c=1.0, kappa=1.0, classes=("C1",), prevalence=0.05),
        np.random.default_rng(5),
    )
    ev = {e["onset"]: e["intensity"] for e in toy_bank["motion"][0]["events"]}
    med = np.median(list(ev.values()))
    for a in rec.artifacts:
        assert a["scale"] == pytest.approx(ev[a["event_onset"]] / med)


def test_level_shift_labels_only_the_transition(toy_bank, cfg):
    """C3 leaves a lasting offset, but only the transition windows are labelled."""
    gen = Generator(toy_bank["library"], cfg)
    rec = gen.generate(
        "x",
        toy_bank["motion"][0],
        toy_bank["phys"][1],
        Condition(c=0.0, classes=("C3",), prevalence=0.01),
        np.random.default_rng(6),
    )
    offset = np.abs(rec.y - rec.y0) >= cfg["generator"]["label_threshold"]
    assert offset.mean() > rec.mask.mean()  # the offset persists after the labelled part
    assert rec.labels.mean() < 0.05
