"""This file belongs to the third step, measuring how real signals behave. It splits an
electrodermal signal into its slow level and its fast responses, because the generator
needs to know both separately.

The method is cvxEDA by Greco et al. (2016), a published and freely available method
that is also offered in open toolboxes for electrodermal activity, such as NeuroKit2
(Makowski et al., 2021). Posada-Quintero and Chon (2020) note in their review that
different decomposition methods give different results, so the program uses one
documented method with the constants of its authors. One choice departs from the
authors: the speed of the response's decay is chosen by the best value of the method's
own objective, because a test on synthetic signals with a known answer showed that the
authors' criterion did not find the right value. The outcome is the slow level, the fast
part and the remainder of every recording.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from scipy.signal import lfilter


def bateman(t: np.ndarray, tau0: float, tau1: float = 0.7) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    tp = np.maximum(t, 0.0)
    return np.where(t >= 0, np.exp(-tp / tau0) - np.exp(-tp / tau1), 0.0)


def tPeak(tau0: float, tau1: float = 0.7) -> float:
    return float(np.log(tau0 / tau1) * tau0 * tau1 / (tau0 - tau1))


def hMax(tau0: float, tau1: float = 0.7) -> float:
    return float(bateman(np.array([tPeak(tau0, tau1)]), tau0, tau1)[0])


def armaCoefficients(delta: float, tau0: float, tau1: float) -> np.ndarray:
    a1 = 1.0 / min(tau1, tau0)
    a0 = 1.0 / max(tau1, tau0)
    return np.array(
        [
            (a1 * delta + 2.0) * (a0 * delta + 2.0),
            2.0 * a1 * a0 * delta**2 - 8.0,
            (a1 * delta - 2.0) * (a0 * delta - 2.0),
        ]
    ) / ((a1 - a0) * delta**2)


def splineBasis(n: int, delta: float, knotS: float = 10.0) -> sparse.csc_matrix:
    k = int(round(knotS / delta))
    knots = np.arange(0, n + k // 2, k)
    tri = np.concatenate((np.arange(1.0, k), np.arange(k, 0.0, -1.0)))
    spl = np.convolve(tri, tri, "full")
    spl /= spl.max()
    offs = np.arange(-(len(spl) // 2), (len(spl) + 1) // 2)
    rows, cols, vals = [], [], []
    for c, kn in enumerate(knots):
        idx = offs + kn
        ok = (idx >= 0) & (idx < n)
        rows.append(idx[ok])
        cols.append(np.full(ok.sum(), c))
        vals.append(spl[ok])
    return sparse.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, len(knots))
    )


class Phasic:

    def __init__(self, ar: np.ndarray):
        self.ar = ar
        self.ma = np.array([1.0, 2.0, 1.0])

    def q(self, p: np.ndarray) -> np.ndarray:
        return lfilter([1.0], self.ar, p)

    def __call__(self, p: np.ndarray) -> np.ndarray:
        return lfilter(self.ma, [1.0], self.q(p))

    def T(self, r: np.ndarray) -> np.ndarray:  # noqa: N802
        z = lfilter(self.ma, [1.0], r[::-1])
        return lfilter([1.0], self.ar, z)[::-1]


def numpySolve(y, delta, tau0, tau1, alpha, gamma, knotS, maxIter, tol, p0=None):
    n = y.size
    K = Phasic(armaCoefficients(delta, tau0, tau1))
    B = splineBasis(n, delta, knotS)
    C = np.column_stack([np.ones(n), np.arange(1.0, n + 1.0) / n])
    W = sparse.hstack([B, sparse.csc_matrix(C)]).tocsc()
    nB = B.shape[1]
    G = (W.T @ W).toarray()
    G[np.arange(nB), np.arange(nB)] += gamma
    chol = cho_factor(G)

    def tonic(res):
        return cho_solve(chol, W.T @ res)

    v = np.random.default_rng(0).standard_normal(n)
    for iteration in range(50):
        v = K.T(K(v))
        v /= np.linalg.norm(v)
    L = float(np.dot(v, K.T(K(v)))) * 1.05

    def obj(p):
        z = tonic(y - K(p))
        e = y - K(p) - W @ z
        return 0.5 * e @ e + alpha * p.sum() + 0.5 * gamma * z[:nB] @ z[:nB], z

    x = np.zeros(n) if p0 is None else np.maximum(p0, 0.0)
    xPrev, yPt, tK, fPrev = x.copy(), x.copy(), 1.0, np.inf
    it = 0
    for it in range(maxIter):
        res = y - K(yPt)
        z = tonic(res)
        grad = -K.T(res - W @ z)
        x = np.maximum(yPt - (grad + alpha) / L, 0.0)
        if np.dot(yPt - x, x - xPrev) > 0:
            tK = 1.0
        tNext = 0.5 * (1 + np.sqrt(1 + 4 * tK**2))
        yPt = x + ((tK - 1) / tNext) * (x - xPrev)
        xPrev, tK = x, tNext
        if it % 50 == 49:
            fNew, unused = obj(x)
            if abs(fPrev - fNew) <= tol * max(1.0, abs(fNew)):
                break
            fPrev = fNew
    p = x
    f, z = obj(p)
    r = K(p)
    t = W @ z
    g = -K.T(y - r - t) + alpha
    kkt = float(np.max(np.abs(np.where(p > 0, g, np.minimum(g, 0.0)))))
    return {"r": r, "p": p, "t": t, "e": y - r - t, "obj": f, "iterations": it + 1, "kkt": kkt}


def referenceModule():
    try:
        import cvxEDA as mod  # type: ignore  # noqa: N813

        return mod
    except ImportError:
        import cvxeda as mod  # type: ignore

        return mod


def decompose(
    y: np.ndarray,
    delta: float,
    tau0: float,
    tau1: float = 0.7,
    alpha: float = 8e-4,
    gamma: float = 1e-2,
    knotS: float = 10.0,
    backend: str = "auto",
    maxIter: int = 20000,
    tol: float = 1e-10,
    p0: np.ndarray | None = None,
) -> dict:
    y = np.asarray(y, dtype=float)
    if backend in ("auto", "cvxeda"):
        try:
            cvxeda = referenceModule()
            r, p, t, unusedL, unusedD, e, obj = cvxeda.cvxEDA(
                y, delta, tau0=tau0, tau1=tau1, delta_knot=knotS, alpha=alpha, gamma=gamma
            )
            return {"r": r, "p": p, "t": t, "e": e, "obj": float(np.ravel(obj)[0]), "backend": "cvxeda"}
        except ImportError:
            if backend == "cvxeda":
                raise
    out = numpySolve(y, delta, tau0, tau1, alpha, gamma, knotS, maxIter, tol, p0)
    out["backend"] = "numpy"
    return out
