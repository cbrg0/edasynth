"""This file is the fourth step and the core of the study. It builds synthetic wristband
recordings in which the cause of every artifact is known.

Real recordings cannot answer the question of the study, because nobody recorded which
artifacts movement caused. The generator therefore assembles new recordings from real
material. Each synthetic person takes the real movements of one participant and the
measured physiology of another, so that any link between movement and skin signal exists
only where the generator puts it. The artifacts copy the shapes of real episodes marked
by experts. Several values come from the literature: responses that follow a stimulus
start one to five seconds later (Boucsein et al., 2012, as summarised by Posada-Quintero
& Chon, 2020), the signal is limited to the valid range of Kleckner et al. (2018), and
3.8 % of windows contain an artifact, as in the labels of Zhang et al. (2017). Applying
the delay of a response to movement is an assumption, because the source describes
responses to stimuli and not to movement. The strength of the link between movement and
artifacts is not fixed at all: it is varied over the scenarios, because it is unknown.

The outcome is a synthetic recording with acceleration, electrodermal activity and a
true label for every window.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import fftconvolve

from .data import WIN, stageType
from .decompose import bateman, hMax, tPeak

UP = 8


@dataclass
class Condition:

    c: float = 0.5
    kappa: float = 1.0
    deltaMax: float = 1.0
    piR: float = 0.0
    classes: tuple = ("C1", "C2", "C3", "C5")
    classWeights: tuple | None = None
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

    def __init__(self, library: dict, cfg: dict):
        self.lib = library
        self.cfg = cfg
        self.g = cfg["generator"]
        mode = self.g["c4"]["gain"]
        if mode not in ("control", "measured"):
            raise ValueError(f"unknown C4 gain mode {mode}")
        self.c4Gain = library.get("Gcontrol" if mode == "control" else "G", float("nan"))
        self.responseJitter = None

    def clean(self, motion: dict, phys: dict, rng: np.random.Generator, piR: float):
        fs = motion["fs"]
        n = motion["acc"].shape[0]
        nf = n * UP
        fsF = fs * UP
        stage = motion["stage"]
        tau0, tau1 = phys["tau0"], self.g["tau1"]
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
        s = np.maximum(s, self.cfg["estimation"]["episodes"]["levelFloor"])
        imp = np.zeros(nf)
        responses = []
        areas = np.asarray(phys["areas"], float)
        for a, b in segments(stage):
            typ = stageType(str(stage[a]))
            lam = toNumber(phys["rates"].get(typ, {}).get("ratePerS"))
            if not np.isfinite(lam):
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
        else:
            j = self.responseJitter
            rF = np.zeros(nf)
            for r in responses:
                t0k = tau0 * (1 + rng.uniform(-j["tau0Rel"], j["tau0Rel"]))
                t1k = rng.uniform(*j["tau1Range"])
                hk = bateman(np.arange(0, self.g["kernelTauMultiple"] * t0k, 1 / fsF), t0k, t1k)
                hk = hk / hMax(t0k, t1k) * hm
                idx = int(round(r["t"] * fsF))
                if idx < nf:
                    seg = hk[: nf - idx] * r["area"]
                    rF[idx : idx + seg.size] += seg
        clean8 = blockMean(s + rF)
        y0 = clean8 + self.noise(phys, clean8.size, rng)
        return y0, responses

    def noise(self, phys: dict, n: int, rng: np.random.Generator) -> np.ndarray:
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
        i0 = min(int(t0 * fs), y0.size - 1)
        w = int(round(self.g["levelWindowS"] * fs))
        lvl = float(np.median(y0[max(0, i0 - w) : i0 + 1]))
        return max(lvl, self.cfg["estimation"]["episodes"]["levelFloor"])

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
        tr = np.asarray(prm["trace"], float)
        i0 = min(int(round(t0 * fs)), n - 1)
        i1 = min(n, i0 + tr.size)
        rel = np.zeros(n)
        rel[i0:i1] = tr[: i1 - i0]
        rel = np.maximum(rel * scale, -self.g["rhoMax"])
        piece = y0 * rel if variant == "V2" else lvl * rel
        sup = np.zeros(n, bool)
        sup[i0:i1] = np.abs(piece[i0:i1]) >= self.g["labelThreshold"]
        return piece, sup

    def c3Event(self, prm, t0, scale, n, fs, lvl):
        i0 = min(int(round(t0 * fs)), n - 1)
        if prm.get("trace") is not None:
            course = np.asarray(prm["trace"], float)
        else:
            k = max(1, int(np.ceil(prm["w"] * fs)))
            course = prm["relDelta"] * np.arange(1, k + 1) / k
        course = course[: n - i0] * scale
        new = prm["relDelta"] * scale
        return {"i0": i0, "course": course, "new": new, "lvl": lvl}

    def contactCourse(self, c3, n, y0, variant):
        lo = -self.g["rhoMax"]
        sup = np.zeros(n, bool)
        state = np.zeros(n)
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
        support = np.zeros(n, bool)
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

        c3 = []
        c3sup = np.zeros(n, bool)

        def positiveWindows():
            m = support | c3sup | (replace & (np.abs(replaceVal - y0) >= thr))
            return m, windowLabels(m, nw, cond.labelRule)

        mask, lab = positiveWindows()
        tries = 0
        while lab.mean() < cond.prevalence and tries < self.g["maxArtifactsPerParticipant"]:
            tries += 1
            cls = cond.classes[rng.choice(len(cond.classes), p=shares)]
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
    return float("nan") if x is None else float(x)


def windowLabels(mask: np.ndarray, nw: int, rule: str) -> np.ndarray:
    m = mask[: nw * WIN].reshape(nw, WIN)
    if rule == "any":
        return m.any(1).astype(int)
    if rule == "majority":
        return (m.sum(1) > WIN / 2).astype(int)
    raise ValueError(rule)


def fastAcceleration(acc: np.ndarray, fs: float, smoothS: float) -> np.ndarray:
    k = int(round(smoothS * fs))
    mag = np.sqrt((acc**2).sum(1))
    return mag - np.convolve(mag, np.ones(k) / k, "same")


def c4ControlGain(motion: list[dict], smoothS: float, targetUs: float) -> float:
    per = []
    for m in motion:
        fast = np.abs(fastAcceleration(m["acc"], m["fs"], smoothS))
        per += [float(np.median(fast[e["onset"] : e["end"]])) for e in m["events"] if e["end"] > e["onset"]]
    med = float(np.median(per)) if per else 0.0
    return targetUs / med if med > 0 else float("nan")


def assignDonors(motion: list[dict], phys: list[dict], rng: np.random.Generator, maxTries: int = 10000) -> np.ndarray:
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
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p
