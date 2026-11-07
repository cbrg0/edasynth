"""This file splits an electrodermal signal into its slow level and its fast responses,
which is needed to measure how real signals behave.

It uses the published method cvxEDA (Greco et al., 2016). The file estimation.py calls
it for every real recording, and the parts it returns tell the generator how large and
how frequent the responses are and how the slow level moves.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from scipy.signal import lfilter


def bateman(t: np.ndarray, tau0: float, tau1: float = 0.7) -> np.ndarray:
    """Return the Bateman function h(t) = (exp(-t/tau0) - exp(-t/tau1)) for t >= 0 and zero
    before, as in equation 4 of Greco et al. (2016).
    """
    t = np.asarray(t, dtype=float)
    tp = np.maximum(t, 0.0)  # This prevents exp(-t/tau) from overflowing for negative t.
    return np.where(t >= 0, np.exp(-tp / tau0) - np.exp(-tp / tau1), 0.0)


def tPeak(tau0: float, tau1: float = 0.7) -> float:
    """Return the time at which the Bateman function reaches its maximum, ln(tau0/tau1)
    tau0 tau1 / (tau0 - tau1).
    """
    return float(np.log(tau0 / tau1) * tau0 * tau1 / (tau0 - tau1))


def hMax(tau0: float, tau1: float = 0.7) -> float:
    """Return the height of the Bateman function at its maximum."""
    return float(bateman(np.array([tPeak(tau0, tau1)]), tau0, tau1)[0])


def armaCoefficients(delta: float, tau0: float, tau1: float) -> np.ndarray:
    """Return the denominator of the Bateman filter after the bilinear transform, as in
    equations 5 to 7 of the paper.
    """
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
    """Return the columns of a smooth cubic curve with joints every knotS seconds.

    Each cubic piece is obtained by convolving two triangles whose half-width equals the
    distance between joints. The columns are centred on the joints 0, k, 2k and so on.
    """
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
    """Apply K p = M A^-1 p and its transpose by recursive filtering."""

    def __init__(self, ar: np.ndarray):
        self.ar = ar
        self.ma = np.array([1.0, 2.0, 1.0])

    def q(self, p: np.ndarray) -> np.ndarray:
        return lfilter([1.0], self.ar, p)

    def __call__(self, p: np.ndarray) -> np.ndarray:
        return lfilter(self.ma, [1.0], self.q(p))

    def T(self, r: np.ndarray) -> np.ndarray:  # noqa: N802 - the name stands for the transposed operator
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

    def tonic(res):  # the exact best (l, d) for the remainder y - K p
        return cho_solve(chol, W.T @ res)

    # The Lipschitz constant of p -> K^T K p is found by power iteration.
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
        # The acceleration is restarted when the gradient points back (O'Donoghue &
        # Candes, 2015).
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
    # Optimality for p >= 0: the gradient is about 0 where p > 0 and not negative where
    # p = 0.
    g = -K.T(y - r - t) + alpha
    kkt = float(np.max(np.abs(np.where(p > 0, g, np.minimum(g, 0.0)))))
    return {"r": r, "p": p, "t": t, "e": y - r - t, "obj": f, "iterations": it + 1, "kkt": kkt}


def referenceModule():
    """Return the authors' cvxEDA code if the user has put it on the Python path, and None
    otherwise.
    """
    try:
        import cvxEDA as mod  # type: ignore  # noqa: N813 - this is the name of the authors' file

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
    """Decompose a signal that is already scaled to mean 0 and standard deviation 1, and
    return r, p, t, e, the value of the objective and the solver used.
    """
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
