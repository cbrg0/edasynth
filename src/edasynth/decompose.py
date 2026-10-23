"""Bateman response and the cvxEDA decomposition (Greco et al., 2016).

Model: y = r + t + e, with phasic r = M A^-1 p (ARMA form of the Bateman
response driven by a sparse non-negative input p), tonic t = B l + C d (cubic
B-spline with 10-s knot spacing plus offset and linear trend) and residual e.
The estimate minimises

    0.5 * ||M A^-1 p + B l + C d - y||^2 + alpha * sum(p) + 0.5 * gamma * ||l||^2,  p >= 0.

Two backends solve the same problem:

``cvxeda``  the authors' reference implementation, file ``cvxEDA.py`` from
            https://github.com/lciti/cvxEDA (needs cvxopt). It is not copied
            here: put the file on the Python path yourself (see README) and
            check its licence there. Imported as ``cvxEDA`` (or ``cvxeda``).
``numpy``   a projected accelerated-gradient solver written for this project,
            dependency-free. ``tests/test_decompose.py`` checks that it recovers
            known inputs and, when ``cvxeda`` is installed, that it reaches the same
            objective value as the reference.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from scipy.signal import lfilter


# ----------------------------------------------------------------- Bateman
def bateman(t: np.ndarray, tau0: float, tau1: float = 0.7) -> np.ndarray:
    """h(t) = (exp(-t/tau0) - exp(-t/tau1)) u(t)  (Greco et al., 2016, eq. 4)."""
    t = np.asarray(t, dtype=float)
    tp = np.maximum(t, 0.0)  # avoid overflow of exp(-t/tau) for negative t
    return np.where(t >= 0, np.exp(-tp / tau0) - np.exp(-tp / tau1), 0.0)


def t_peak(tau0: float, tau1: float = 0.7) -> float:
    """Time of the maximum of the Bateman function: ln(tau0/tau1) tau0 tau1 / (tau0 - tau1)."""
    return float(np.log(tau0 / tau1) * tau0 * tau1 / (tau0 - tau1))


def h_max(tau0: float, tau1: float = 0.7) -> float:
    """Height of the Bateman function at its peak."""
    return float(bateman(np.array([t_peak(tau0, tau1)]), tau0, tau1)[0])


# ----------------------------------------------------------------- operators
def arma_coefficients(delta: float, tau0: float, tau1: float) -> np.ndarray:
    """Denominator of the bilinear-transformed Bateman filter (paper eqs. 5-7)."""
    a1 = 1.0 / min(tau1, tau0)
    a0 = 1.0 / max(tau1, tau0)
    return np.array(
        [
            (a1 * delta + 2.0) * (a0 * delta + 2.0),
            2.0 * a1 * a0 * delta**2 - 8.0,
            (a1 * delta - 2.0) * (a0 * delta - 2.0),
        ]
    ) / ((a1 - a0) * delta**2)


def spline_basis(n: int, delta: float, knot_s: float = 10.0) -> sparse.csc_matrix:
    """Cubic B-spline regressors with knots every ``knot_s`` seconds.

    The cubic B-spline is obtained by convolving two triangles of half-width
    equal to the knot spacing; columns are centred on the knots 0, k, 2k, ...
    """
    k = int(round(knot_s / delta))
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


class _Phasic:
    """K p = M A^-1 p and its transpose, applied by recursive filtering."""

    def __init__(self, ar: np.ndarray):
        self.ar = ar
        self.ma = np.array([1.0, 2.0, 1.0])

    def q(self, p: np.ndarray) -> np.ndarray:
        return lfilter([1.0], self.ar, p)

    def __call__(self, p: np.ndarray) -> np.ndarray:
        return lfilter(self.ma, [1.0], self.q(p))

    def T(self, r: np.ndarray) -> np.ndarray:  # noqa: N802 - transpose operator
        z = lfilter(self.ma, [1.0], r[::-1])
        return lfilter([1.0], self.ar, z)[::-1]


def _numpy_solve(y, delta, tau0, tau1, alpha, gamma, knot_s, max_iter, tol, p0=None):
    n = y.size
    K = _Phasic(arma_coefficients(delta, tau0, tau1))
    B = spline_basis(n, delta, knot_s)
    C = np.column_stack([np.ones(n), np.arange(1.0, n + 1.0) / n])
    W = sparse.hstack([B, sparse.csc_matrix(C)]).tocsc()
    nB = B.shape[1]
    G = (W.T @ W).toarray()
    G[np.arange(nB), np.arange(nB)] += gamma
    chol = cho_factor(G)

    def tonic(res):  # exact minimiser over (l, d) given the residual y - K p
        return cho_solve(chol, W.T @ res)

    # Lipschitz constant of p -> K^T K p by power iteration
    v = np.random.default_rng(0).standard_normal(n)
    for _ in range(50):
        v = K.T(K(v))
        v /= np.linalg.norm(v)
    L = float(np.dot(v, K.T(K(v)))) * 1.05

    def obj(p):
        z = tonic(y - K(p))
        e = y - K(p) - W @ z
        return 0.5 * e @ e + alpha * p.sum() + 0.5 * gamma * z[:nB] @ z[:nB], z

    x = np.zeros(n) if p0 is None else np.maximum(p0, 0.0)
    x_prev, y_pt, t_k, f_prev = x.copy(), x.copy(), 1.0, np.inf
    it = 0
    for it in range(max_iter):
        res = y - K(y_pt)
        z = tonic(res)
        grad = -K.T(res - W @ z)
        x = np.maximum(y_pt - (grad + alpha) / L, 0.0)
        # gradient-based adaptive restart (O'Donoghue & Candes, 2015)
        if np.dot(y_pt - x, x - x_prev) > 0:
            t_k = 1.0
        t_next = 0.5 * (1 + np.sqrt(1 + 4 * t_k**2))
        y_pt = x + ((t_k - 1) / t_next) * (x - x_prev)
        x_prev, t_k = x, t_next
        if it % 50 == 49:
            f_new, _ = obj(x)
            if abs(f_prev - f_new) <= tol * max(1.0, abs(f_new)):
                break
            f_prev = f_new
    p = x
    f, z = obj(p)
    r = K(p)
    t = W @ z
    # optimality check (KKT for p >= 0): gradient ~0 where p > 0, gradient >= 0 where p = 0
    g = -K.T(y - r - t) + alpha
    kkt = float(np.max(np.abs(np.where(p > 0, g, np.minimum(g, 0.0)))))
    return {"r": r, "p": p, "t": t, "e": y - r - t, "obj": f, "iterations": it + 1, "kkt": kkt}


def _reference_module():
    """The reference cvxEDA implementation, if the user has put it on the path."""
    try:
        import cvxEDA as mod  # type: ignore  # noqa: N813 - name of the authors' file

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
    knot_s: float = 10.0,
    backend: str = "auto",
    max_iter: int = 20000,
    tol: float = 1e-10,
    p0: np.ndarray | None = None,
) -> dict:
    """Decompose an (already z-scored) signal; returns r, p, t, e, obj and the backend used."""
    y = np.asarray(y, dtype=float)
    if backend in ("auto", "cvxeda"):
        try:
            cvxeda = _reference_module()
            r, p, t, _l, _d, e, obj = cvxeda.cvxEDA(
                y, delta, tau0=tau0, tau1=tau1, delta_knot=knot_s, alpha=alpha, gamma=gamma
            )
            return {"r": r, "p": p, "t": t, "e": e, "obj": float(np.ravel(obj)[0]), "backend": "cvxeda"}
        except ImportError:
            if backend == "cvxeda":
                raise
    out = _numpy_solve(y, delta, tau0, tau1, alpha, gamma, knot_s, max_iter, tol, p0)
    out["backend"] = "numpy"
    return out
