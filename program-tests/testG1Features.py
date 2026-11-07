"""G1: features are computed as documented, causally, and reproduce the published ones."""

import numpy as np
import pytest
from conftest import ZHANG, needsZhang
from scipy import stats

from edasynth.features import derivatives, f0, f0Names, fa, fa4, faNames, fb, haarDetails, windows


def testWindowHas40Samples():
    assert windows(np.arange(100.0)).shape == (2, 40)


def testFeatureNameCounts():
    assert len(f0Names) == 24 and len(faNames) == 96


def testHaarAgainstHandComputation():
    x = np.arange(8.0)
    c1, c2, c3 = haarDetails(np.tile(x, 5)[:40])
    assert np.allclose(c1[:4], (x[0::2] - x[1::2]) / np.sqrt(2))
    assert c3.shape == (5,)


def testHaarMatchesPywaveletsIfInstalled():
    pywt = pytest.importorskip("pywt")
    w = np.random.default_rng(0).standard_normal(40)
    ref0, ref3, ref2, ref1 = pywt.wavedec(w, "haar", level=3)
    c1, c2, c3 = haarDetails(w)
    assert np.allclose(c1, ref1) and np.allclose(c2, ref2) and np.allclose(c3, ref3)


def testCentralDerivativeEqualsZhangFormula():
    w = np.random.default_rng(1).standard_normal(40)
    d1, d2 = derivatives(w, "central")
    zhang = (w[1:-1] + w[2:]) / 2.0 - (w[1:-1] + w[:-2]) / 2.0
    assert np.allclose(d1, zhang)
    assert np.allclose(d2, w[2:] - 2 * w[1:-1] + w[:-2])


def testFa4IsFirstFourOfFa():
    acc = np.random.default_rng(2).standard_normal((400, 3))
    assert np.allclose(fa(acc)[:, :4], fa4(acc))


def testFbIsCausal():
    """Changing samples after the end of window w leaves the features of w unchanged."""
    rng = np.random.default_rng(3)
    n = 8 * 300
    y = 2 + np.cumsum(rng.standard_normal(n)) * 0.001
    acc = rng.standard_normal((n, 3)) * 0.1
    a = fb(y, acc, 8.0, 0.001, 0.01)
    w = 20
    y2, acc2 = y.copy(), acc.copy()
    y2[(w + 1) * 40 :] += rng.standard_normal(n - (w + 1) * 40)
    acc2[(w + 1) * 40 :] += 5
    b = fb(y2, acc2, 8.0, 0.001, 0.01)
    assert np.allclose(a[: w + 1], b[: w + 1])
    assert not np.allclose(a[w + 1 :], b[w + 1 :])


def testB3IsFiniteOnFlatSignal():
    n = 8 * 120
    out = fb(np.full(n, 3.0), np.tile([0, 0, 1.0], (n, 1)), 8.0, 0.01, 0.01)
    assert np.all(np.isfinite(out))


def testB1FindsPlantedLag():
    rng = np.random.default_rng(4)
    n = 8 * 120
    acc = np.tile([0, 0, 1.0], (n, 1)) + 0.001 * rng.standard_normal((n, 3))
    y = np.full(n, 2.0) + 0.0001 * rng.standard_normal(n)
    spikes = [330, 341, 349, 362, 371]  # irregular: a periodic train would alias the lag
    for s in spikes:
        acc[s] += [0.5, 0, 0]
        y[s + 16 :] += 0.05  # step 2 s after each movement
    out = fb(y, acc, 8.0, 1e-4, 1e-3)
    w = 9  # window containing samples 360-399
    assert out[w, 0] > 0.5
    assert out[w, 1] == pytest.approx(2.0, abs=0.25)


@needsZhang
def testAmplitudeAndSecondDerivativeReproducePublishedAwwFeatures():
    """Recomputed AWW walking features vs the published matrix (z-scored per recording)."""
    from edasynth.data import loadAww

    recs = loadAww(ZHANG)[:5]
    ours = np.vstack([stats.zscore(f0(r.eda, "forward"), axis=0) for r in recs])
    ref = np.loadtxt(ZHANG / "AlanWalksWales" / "AWW_walk_eda.csv", delimiter=",")
    for k in list(range(0, 12)):
        r = np.corrcoef(ours[:, k], ref[:, k])[0, 1]
        assert r > 0.999, (f0Names[k], r)


@needsZhang
def testPublishedFirstDerivativeIsForwardNotCentral():
    from edasynth.data import loadAww

    recs = loadAww(ZHANG)[:5]
    ref = np.loadtxt(ZHANG / "AlanWalksWales" / "AWW_walk_eda.csv", delimiter=",")
    cen = np.vstack([stats.zscore(f0(r.eda, "central"), axis=0) for r in recs])
    rs = [np.corrcoef(cen[:, k], ref[:, k])[0, 1] for k in range(4, 8)]
    assert min(rs) < 0.995
