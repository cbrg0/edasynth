"""Semi-synthetic generator (docs/METHOD.md section 4, including the movement-artifact link).

Each synthetic participant combines
  * the real accelerometer record, protocol stages and motion events of a
    *motion donor* u(j), and
  * the physiology (tau0, response rates and sizes, tonic shape, noise,
    quantisation) of a different *physiology donor* d(j) != u(j), chosen among
    the donors whose EDA level is valid.
Because the donors are different people, every link between movement and the
EDA signal in the output is created by the generator's parameters and is
therefore known exactly.

Clean signal   y0 = s + r + eps            (s, r on a 64-Hz grid, averaged to 8 Hz; eps =
                                           blocks of the donor's decomposition residual)
Artifact size  relative template size x local level of y0 before the artifact
Artifacts      y  = y0 + sum of all artifacts (variant V1, additive)
               y  = y0 * (1 - rho)          (variant V2, classes C2/C3 only)
True label     sample corrupted if the artifact moves it by >= threshold (for the
               level shift C3 only the transition counts, not the lasting offset);
               window positive if any (or, as a variant, most) samples are corrupted.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import fftconvolve

from .data import WIN, stageType
from .decompose import bateman, hMax, tPeak

UP = 8  # fine grid = 8 x 8 Hz = 64 Hz


@dataclass
class Condition:
    """One cell of the grid (docs/METHOD.md section 8)."""

    c: float = 0.5  # share of artifacts caused by movement
    kappa: float = 1.0  # amplitude ~ (event intensity / median intensity)^kappa
    deltaMax: float = 1.0  # maximum delay after motion onset, s
    piR: float = 0.0  # share of motion events followed by a real response
    classes: tuple = ("C1", "C2", "C3", "C5")  # composition
    classWeights: tuple | None = None  # None = measured UTD shares
    variant: str = "V1"
    prevalence: float = 0.038
    labelRule: str = "any"

    def key(self) -> str:
        cl = "+".join(self.classes)
        return (
            f"c={self.c}|kappa={self.kappa}|dmax={self.deltaMax}|piR={self.piR}|"
            f"cls={cl}|{self.variant}|prev={self.prevalence}|{self.labelRule}"
        )


@dataclass
class SynthRecord:
    pid: str
    fs: float
    y: np.ndarray
    y0: np.ndarray
    acc: np.ndarray
    stage: np.ndarray
    mask: np.ndarray
    labels: np.ndarray
    artifacts: list = field(default_factory=list)
    responses: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)


def segments(stage: np.ndarray):
    starts = np.r_[0, np.where(stage[1:] != stage[:-1])[0] + 1]
    ends = np.r_[starts[1:], stage.size]
    return list(zip(starts, ends, strict=True))


def blockMean(x: np.ndarray, k: int = UP) -> np.ndarray:
    n = x.size // k
    return x[: n * k].reshape(n, k).mean(1)


def trapezoid(nf: int, fsF: float, t0: float, rise: float, plateau: float, fall: float, height: float):
    t = np.arange(nf) / fsF - t0
    a = np.zeros(nf)
    r = max(rise, 1 / fsF)
    f = max(fall, 1 / fsF)
    up = (t >= 0) & (t < r)
    a[up] = t[up] / r
    pl = (t >= r) & (t < r + plateau)
    a[pl] = 1.0
    dn = (t >= r + plateau) & (t < r + plateau + f)
    a[dn] = 1 - (t[dn] - r - plateau) / f
    return height * a


def stepShape(nf: int, fsF: float, t0: float, w: float, height: float):
    t = np.arange(nf) / fsF - t0
    w = max(w, 1 / fsF)
    a = np.clip(t / w, 0, 1)
    return height * a


class Generator:
    """Builds synthetic participants from estimated parameters."""

    def __init__(self, library: dict, cfg: dict):
        self.lib = library  # {"C1": [params...], ..., "classShares": {...}, "G": float}
        self.cfg = cfg
        self.g = cfg["generator"]

    # ------------------------------------------------------------ clean part
    def clean(self, motion: dict, phys: dict, rng: np.random.Generator, piR: float):
        fs = motion["fs"]
        n = motion["acc"].shape[0]
        nf = n * UP
        fsF = fs * UP
        stage = motion["stage"]
        tau0, tau1 = phys["tau0"], self.g["tau1"]
        # tonic: one fitted quadruple per protocol stage, stitched without jumps
        s = np.zeros(nf)
        level = None
        for a, b in segments(stage):
            typ = stageType(str(stage[a]))
            pool = [t for t in phys["tonic"] if t["stageType"] == typ] or phys["tonic"]
            q4 = pool[rng.integers(len(pool))]
            tau = np.arange((b - a) * UP) / fsF
            seg = q4["b"] * tau + q4["AS"] * np.sin(2 * np.pi * tau / q4["TT"] + rng.uniform(0, 2 * np.pi))
            start = q4["s0"] if level is None else level
            seg = seg - seg[0] + start
            s[a * UP : b * UP] = seg
            level = seg[-1]
        # responses: spontaneous (Poisson) and motion-evoked (share piR)
        imp = np.zeros(nf)
        responses = []
        areas = np.asarray(phys["areas"], float)
        for a, b in segments(stage):
            typ = stageType(str(stage[a]))
            lam = toNumber(phys["rates"].get(typ, {}).get("ratePerS"))
            if not np.isfinite(lam):  # donor has no clean time of this stage type: use its other rate
                other = [
                    toNumber(v["ratePerS"]) for v in phys["rates"].values() if np.isfinite(toNumber(v["ratePerS"]))
                ]
                lam = other[0] if other else 0.0
                self.rateFallbacks = getattr(self, "rateFallbacks", 0) + 1
            dur = (b - a) / fs
            k = rng.poisson(lam * dur) if np.isfinite(lam) else 0
            for tK in np.sort(rng.uniform(a / fs, b / fs, k)):
                area = float(areas[rng.integers(areas.size)])
                responses.append({"t": float(tK), "area": area, "origin": "spontaneous"})
        for ev in motion["events"]:
            if rng.random() < piR:
                lat = rng.uniform(*self.g["motionResponseLatencyS"])
                tK = ev["onset"] / fs + lat
                if tK < n / fs:
                    area = float(areas[rng.integers(areas.size)])
                    responses.append({"t": float(tK), "area": area, "origin": "motion"})
        hm, tp = hMax(tau0, tau1), tPeak(tau0, tau1)
        for r in responses:
            idx = int(round(r["t"] * fsF))
            if idx < nf:
                imp[idx] += r["area"]
            r["amp"] = r["area"] * hm
            r["peakSample"] = int((r["t"] + tp) * fs)
        h = bateman(np.arange(0, self.g["kernelTauMultiple"] * tau0, 1 / fsF), tau0, tau1)
        rF = fftconvolve(imp, h)[:nf]
        clean8 = blockMean(s + rF)
        y0 = clean8 + self.noise(phys, clean8.size, rng)
        return y0, responses

    def noise(self, phys: dict, n: int, rng: np.random.Generator) -> np.ndarray:
        """Measurement noise of the clean signal (docs/METHOD.md section 4.2).

        ``residualBootstrap``: blocks of the physiology donor's decomposition
        residual, drawn with replacement and joined; ``gaussian``: independent
        normal noise with the donor's sigma.
        """
        mode = self.g["noise"]
        if mode == "gaussian":
            return rng.normal(0.0, phys["sigma"], n)
        if mode != "residualBootstrap":
            raise ValueError(f"unknown noise model {mode}")
        blocks = np.asarray(phys.get("residual", np.zeros((0, 1))), float)
        if blocks.ndim != 2 or blocks.shape[0] == 0:
            raise ValueError(f"{phys['pid']}: no residual blocks in the bank for the noise model")
        L = blocks.shape[1]
        idx = rng.integers(blocks.shape[0], size=-(-n // L))
        return blocks[idx].ravel()[:n]

    def localLevel(self, y0: np.ndarray, t0: float, fs: float) -> float:
        """Median of the clean signal over the ``levelWindowS`` before t0, floored."""
        i0 = min(int(t0 * fs), y0.size - 1)
        w = int(round(self.g["levelWindowS"] * fs))
        lvl = float(np.median(y0[max(0, i0 - w) : i0 + 1]))
        return max(lvl, self.cfg["estimation"]["episodes"]["levelFloor"])

    # ------------------------------------------------------------ artifacts
    def drawTemplate(self, cls: str, rng: np.random.Generator):
        pool = self.lib.get(cls) or []
        if pool:
            return dict(pool[rng.integers(len(pool))]), "measured"
        fb = self.g["fallbackTemplates"].get(cls)
        if fb is None or not self.g["allowFallbackTemplates"]:
            raise ValueError(f"no measured {cls} episodes and fallback templates disabled")
        out = {}
        for k, v in fb.items():
            out[k] = rng.uniform(*v) if isinstance(v, list) else v
        return out, "fallback"

    def artifactSignal(self, cls, prm, t0, scale, n, fs, y0, acc, ev, variant):
        nf, fsF = n * UP, fs * UP
        # template sizes are relative to the level before the artifact (section 4.3)
        lvl = self.localLevel(y0, t0, fs) if cls in ("C1", "C2", "C3") else 1.0
        if cls == "C1":
            a = blockMean(trapezoid(nf, fsF, t0, prm["tauR"], 0.0, prm["tauF"], prm["relA"] * lvl * scale))
            return a, np.abs(a) >= self.g["labelThreshold"]
        if cls in ("C2", "C3"):
            if cls == "C2":
                shape = trapezoid(nf, fsF, t0, prm["w"], prm["D"], prm["w"], 1.0)
                size = -prm["relDepth"] * lvl * scale
            else:
                shape = stepShape(nf, fsF, t0, prm["w"], 1.0)
                size = prm["relDelta"] * lvl * scale
            shape8 = blockMean(shape)
            if variant == "V2":
                rho = np.clip(abs(size) / lvl, 0.0, self.g["rhoMax"])
                piece = -y0 * rho * shape8
            else:
                piece = size * shape8
            if cls == "C3":
                # The lasting offset after a contact change is not labelled; only the
                # transition is (plus the first sample at the new level), as experts do.
                moving = (shape8 > 0) & (shape8 < 1)
                first = np.argmax(shape8 >= 1) if np.any(shape8 >= 1) else None
                if first is not None:
                    moving[first] = True
                return piece, moving & (abs(size) >= self.g["labelThreshold"])
            return piece, np.abs(piece) >= self.g["labelThreshold"]
        if cls == "C4":
            k = int(round(self.cfg["estimation"]["c4"]["smoothS"] * fs))
            mag = np.sqrt((acc**2).sum(1))
            fast = mag - np.convolve(mag, np.ones(k) / k, "same")
            out = np.zeros(n)
            out[ev["onset"] : ev["end"]] = self.lib["G"] * scale * fast[ev["onset"] : ev["end"]]
            return out, np.abs(out) >= self.g["labelThreshold"]
        raise ValueError(cls)

    def generate(self, pid: str, motion: dict, phys: dict, cond: Condition, rng: np.random.Generator) -> SynthRecord:
        fs = motion["fs"]
        acc = motion["acc"]
        n = acc.shape[0]
        y0, responses = self.clean(motion, phys, rng, cond.piR)
        thr = self.g["labelThreshold"]
        art = np.zeros(n)
        support = np.zeros(n, bool)  # samples labelled as corrupted
        artByCls: dict[str, np.ndarray] = {}
        replace = np.zeros(n, bool)
        replaceVal = np.zeros(n)
        artifacts = []
        events = motion["events"]
        medI = float(np.median([e["intensity"] for e in events])) if events else 1.0
        nw = n // WIN
        shares = cond.classWeights or tuple(self.lib["classShares"].get(c, 0.0) for c in cond.classes)
        shares = np.asarray(shares, float)
        shares = shares / shares.sum() if shares.sum() > 0 else np.full(len(cond.classes), 1 / len(cond.classes))

        def positiveWindows():
            m = support | (replace & (np.abs(replaceVal - y0) >= thr))
            return m, windowLabels(m, nw, cond.labelRule)

        mask, lab = positiveWindows()
        tries = 0
        while lab.mean() < cond.prevalence and tries < self.g["maxArtifactsPerParticipant"]:
            tries += 1
            cls = cond.classes[rng.choice(len(cond.classes), p=shares)]
            # C4 has no template: its shape is G * (fast part of |a|) inside a motion event
            prm, src = ({}, "modelG") if cls == "C4" else self.drawTemplate(cls, rng)
            ev = None
            if cls == "C4":
                if not events:
                    break
                ev = events[rng.integers(len(events))]
                t0, caused = ev["onset"] / fs, True
                scale = (ev["intensity"] / medI) ** cond.kappa
            elif events and rng.random() < cond.c:
                ev = events[rng.integers(len(events))]
                t0 = ev["onset"] / fs + rng.uniform(0, cond.deltaMax)
                caused = True
                scale = (ev["intensity"] / medI) ** cond.kappa
            else:
                t0, caused, scale = rng.uniform(0, n / fs), False, 1.0
            if cls == "C5":
                a = int(t0 * fs)
                b = min(n, a + max(1, int(round(prm["D"] * fs))))
                replace[a:b] = True
                replaceVal[a:b] = 0.0 if prm.get("kind", "zero") == "zero" else self.g["saturationLevel"]
            else:
                piece, sup = self.artifactSignal(cls, prm, t0, scale, n, fs, y0, acc, ev, cond.variant)
                art += piece
                support |= sup
                artByCls[cls] = artByCls.get(cls, np.zeros(n, bool)) | sup
            artifacts.append(
                {
                    "cls": cls,
                    "t0": float(t0),
                    "causedByMotion": caused,
                    "scale": float(scale),
                    "eventOnset": None if ev is None else int(ev["onset"]),
                    "template": src,
                    **prm,
                }
            )
            mask, lab = positiveWindows()
        y = y0 + art
        y[replace] = replaceVal[replace]
        classWindows = {c: windowLabels(a, nw, "any") for c, a in artByCls.items()}
        if replace.any():
            classWindows["C5"] = windowLabels(replace & (np.abs(replaceVal - y0) >= thr), nw, "any")
        q = phys["q"]
        if q > 0:
            y = q * np.round(y / q)
            y0q = q * np.round(y0 / q)
        else:
            y0q = y0
        return SynthRecord(
            pid=pid,
            fs=fs,
            y=y,
            y0=y0q,
            acc=acc,
            stage=motion["stage"],
            mask=mask,
            labels=lab,
            artifacts=artifacts,
            responses=responses,
            meta={
                "motionDonor": motion["pid"],
                "physiologyDonor": phys["pid"],
                "condition": cond.key(),
                "q": q,
                "qAcc": motion.get("qAcc", 0.0),
                "classWindows": classWindows,
            },
        )


def toNumber(x) -> float:
    """Number from the bank; JSON null (missing estimate) becomes NaN."""
    return float("nan") if x is None else float(x)


def windowLabels(mask: np.ndarray, nw: int, rule: str) -> np.ndarray:
    m = mask[: nw * WIN].reshape(nw, WIN)
    if rule == "any":
        return m.any(1).astype(int)
    if rule == "majority":
        return (m.sum(1) > WIN / 2).astype(int)
    raise ValueError(rule)


def assignDonors(motion: list[dict], phys: list[dict], rng: np.random.Generator, maxTries: int = 10000) -> np.ndarray:
    """Physiology donor index for every movement donor (docs/METHOD.md section 4.1).

    Only valid physiology donors are used and nobody is paired with themself.
    When every donor is valid this is a derangement, so each physiology donor is
    used once. Otherwise the valid donors are dealt out in random order, each at
    most one time more often than any other.
    """
    n = len(motion)
    valid = np.array([i for i, p in enumerate(phys) if p.get("validDonor", True)])
    if valid.size == n:
        return derangement(n, rng)
    if valid.size < 2:
        raise ValueError("fewer than two valid physiology donors")
    mPid = np.array([m["pid"] for m in motion])
    pPid = np.array([p["pid"] for p in phys])
    reps = -(-n // valid.size)
    for attempt in range(maxTries):
        perm = np.concatenate([rng.permutation(valid) for r in range(reps)])[:n]
        perm = rng.permutation(perm)
        if not np.any(pPid[perm] == mPid):
            return perm
    raise RuntimeError("could not pair every movement donor with a different physiology donor")


def derangement(n: int, rng: np.random.Generator) -> np.ndarray:
    """Random permutation with no fixed point (physiology donor != motion donor)."""
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p
