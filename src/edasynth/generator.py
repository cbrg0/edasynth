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
        self.lib = library  # {"C1": [params...], ..., "classShares": {...}, "G": float, "Gcontrol": float}
        self.cfg = cfg
        self.g = cfg["generator"]
        mode = self.g["c4"]["gain"]
        if mode not in ("control", "measured"):
            raise ValueError(f"unknown C4 gain mode {mode}")
        self.c4Gain = library.get("Gcontrol" if mode == "control" else "G", float("nan"))
        # set by the G2 check only: per-response variation of the time constants
        self.responseJitter = None

    # ------------------------------------------------------------ clean part
    def clean(self, motion: dict, phys: dict, rng: np.random.Generator, piR: float):
        fs = motion["fs"]
        n = motion["acc"].shape[0]
        nf = n * UP
        fsF = fs * UP
        stage = motion["stage"]
        tau0, tau1 = phys["tau0"], self.g["tau1"]
        # tonic (docs/METHOD.md section 4.2): for each stage of the movement donor the
        # physiology donor's fit of the same stage (else one of the same type) gives the
        # level reached at the end of the stage and the oscillation; the trend runs from
        # the level at the start of the stage to that end level, so the level follows
        # the donor's own course instead of adding up slopes of unrelated stages
        s = np.zeros(nf)
        level = None
        for a, b in segments(stage):
            name = str(stage[a])
            typ = stageType(name)
            pool = [t for t in phys["tonic"] if t.get("stage") == name]
            pool = pool or [t for t in phys["tonic"] if t["stageType"] == typ] or phys["tonic"]
            q4 = pool[rng.integers(len(pool))]
            tau = np.arange((b - a) * UP) / fsF
            start = q4["s0"] if level is None else level
            end = q4["s0"] + q4["b"] * q4["durationS"]
            trend = start + (end - start) * tau / max(tau[-1], 1 / fsF)
            osc = q4["AS"] * np.sin(2 * np.pi * tau / q4["TT"] + rng.uniform(0, 2 * np.pi))
            seg = trend + osc - osc[0]
            s[a * UP : b * UP] = seg
            level = seg[-1]
        # a line plus a sinusoid can dip below what a skin can show; the tonic is kept at or
        # above the minimum valid level (Kleckner et al., 2018, rule 1) [decision]
        s = np.maximum(s, self.cfg["estimation"]["episodes"]["levelFloor"])
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
        if self.responseJitter is None:
            h = bateman(np.arange(0, self.g["kernelTauMultiple"] * tau0, 1 / fsF), tau0, tau1)
            rF = fftconvolve(imp, h)[:nf]
        else:  # each response with its own time constants, so the shape differs from cvxEDA's
            j = self.responseJitter
            rF = np.zeros(nf)
            for r in responses:
                t0k = tau0 * (1 + rng.uniform(-j["tau0Rel"], j["tau0Rel"]))
                t1k = rng.uniform(*j["tau1Range"])
                hk = bateman(np.arange(0, self.g["kernelTauMultiple"] * t0k, 1 / fsF), t0k, t1k)
                # same height as the unjittered response of this size
                hk = hk / hMax(t0k, t1k) * hm
                idx = int(round(r["t"] * fsF))
                if idx < nf:
                    seg = hk[: nf - idx] * r["area"]
                    rF[idx : idx + seg.size] += seg
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
        if cls in ("C1", "C2") and prm.get("trace") is not None:
            return self.transplant(cls, prm, t0, scale, n, fs, y0, lvl, variant)
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
            fast = fastAcceleration(acc, fs, self.cfg["estimation"]["c4"]["smoothS"])
            out = np.zeros(n)
            out[ev["onset"] : ev["end"]] = self.c4Gain * scale * fast[ev["onset"] : ev["end"]]
            return out, np.abs(out) >= self.g["labelThreshold"]
        raise ValueError(cls)

    def transplant(self, cls, prm, t0, scale, n, fs, y0, lvl, variant):
        """Peak or drop with the measured course of a real episode (docs/METHOD.md section 4.3).

        The relative trace, times the movement scaling and bounded below by
        -rhoMax, is multiplied by the local level (V1) or by the signal itself (V2).
        Level shifts (C3) are built by :meth:`contactCourse`.
        """
        tr = np.asarray(prm["trace"], float)
        i0 = min(int(round(t0 * fs)), n - 1)
        i1 = min(n, i0 + tr.size)
        rel = np.zeros(n)
        rel[i0:i1] = tr[: i1 - i0]
        # an artifact can remove at most rhoMax of the signal, also after movement scaling
        rel = np.maximum(rel * scale, -self.g["rhoMax"])
        piece = y0 * rel if variant == "V2" else lvl * rel
        sup = np.zeros(n, bool)
        sup[i0:i1] = np.abs(piece[i0:i1]) >= self.g["labelThreshold"]
        return piece, sup

    def c3Event(self, prm, t0, scale, n, fs, lvl):
        """One level shift as (start, relative course during the change, new relative level)."""
        i0 = min(int(round(t0 * fs)), n - 1)
        if prm.get("trace") is not None:
            course = np.asarray(prm["trace"], float)
        else:  # simple shape: linear ramp of width w to the shift
            k = max(1, int(np.ceil(prm["w"] * fs)))
            course = prm["relDelta"] * np.arange(1, k + 1) / k
        course = course[: n - i0] * scale
        new = prm["relDelta"] * scale
        return {"i0": i0, "course": course, "new": new, "lvl": lvl}

    def contactCourse(self, c3, n, y0, variant):
        """Offset produced by all level shifts, and the samples labelled as corrupted.

        The level after a shift (the contact state) replaces the previous state instead
        of adding to it, so shifts do not accumulate [decision, 1.4.0]. During a shift
        the signal moves from the previous state to the measured course of the episode
        (linear blend over the episode); only that change is labelled, the lasting
        offset afterwards is not. Relative values are bounded below by -rhoMax.
        """
        lo = -self.g["rhoMax"]
        sup = np.zeros(n, bool)
        state = np.zeros(n)  # in uS for V1, relative to y0 for V2
        for e in sorted(c3, key=lambda e: e["i0"]):
            unit = 1.0 if variant == "V2" else e["lvl"]
            i0 = e["i0"]
            course = np.maximum(e["course"], lo) * unit
            i1 = i0 + course.size
            prev = state[i0]
            v = np.arange(1, course.size + 1) / course.size
            during = prev * (1 - v) + course
            state[i0:i1] = during
            state[i1:] = max(e["new"], lo) * unit
            change = np.abs(during - prev) * (y0[i0:i1] if variant == "V2" else 1.0)
            sup[i0:i1] |= change >= self.g["labelThreshold"]
        return (y0 * state if variant == "V2" else state), sup

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

        c3 = []  # level shifts, combined by contactCourse
        c3sup = np.zeros(n, bool)

        def positiveWindows():
            m = support | c3sup | (replace & (np.abs(replaceVal - y0) >= thr))
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
            elif cls == "C3":
                c3.append(self.c3Event(prm, t0, scale, n, fs, self.localLevel(y0, t0, fs)))
                unused, c3sup = self.contactCourse(c3, n, y0, cond.variant)
            else:
                # sizes refer to the level of the clean signal, so artifacts add instead of compounding
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
                    **{k: v for k, v in prm.items() if k != "trace"},
                }
            )
            mask, lab = positiveWindows()
        if c3:
            off, c3sup = self.contactCourse(c3, n, y0, cond.variant)
            art = art + off
            artByCls["C3"] = c3sup
        y = y0 + art
        # together the artifacts can remove at most rhoMax of the signal, and the signal
        # cannot exceed the upper range limit of the sensor (Kleckner et al., 2018, rule 1)
        y = np.clip(y, (1 - self.g["rhoMax"]) * y0, self.g["saturationLevel"])
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


def fastAcceleration(acc: np.ndarray, fs: float, smoothS: float) -> np.ndarray:
    """|a| minus its moving mean over ``smoothS`` seconds (the shape of class C4)."""
    k = int(round(smoothS * fs))
    mag = np.sqrt((acc**2).sum(1))
    return mag - np.convolve(mag, np.ones(k) / k, "same")


def c4ControlGain(motion: list[dict], smoothS: float, targetUs: float) -> float:
    """Gain that gives the median movement event a C4 artifact of median size ``targetUs``.

    For every movement event of every movement donor the median of |fast
    acceleration| inside the event is taken; the gain is the target divided by the
    median of these values (docs/METHOD.md section 4.3). Used because the gain
    measured on UTD is about 0, so C4 would produce no artifact [decision].
    """
    per = []
    for m in motion:
        fast = np.abs(fastAcceleration(m["acc"], m["fs"], smoothS))
        per += [float(np.median(fast[e["onset"] : e["end"]])) for e in m["events"] if e["end"] > e["onset"]]
    med = float(np.median(per)) if per else 0.0
    return targetUs / med if med > 0 else float("nan")


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
