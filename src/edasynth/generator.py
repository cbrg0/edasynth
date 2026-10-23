"""Semi-synthetic generator (docs/METHOD.md section 4, including the movement-artifact link).

Each synthetic participant combines
  * the real accelerometer record, protocol stages and motion events of a
    *motion donor* u(j), and
  * the physiology (tau0, response rates and sizes, tonic shape, noise,
    quantisation) of a different *physiology donor* d(j) != u(j).
Because the donors are different people, every link between movement and the
EDA signal in the output is created by the generator's parameters and is
therefore known exactly.

Clean signal   y0 = s + r + eps            (built on a 64-Hz grid, averaged to 8 Hz)
Artifacts      y  = y0 + sum_m A_m          (variant V1, additive)
               y  = y0 * (1 - rho)          (variant V2, classes C2/C3 only)
True label     sample corrupted if the artifact moves it by >= threshold (for the
               level shift C3 only the transition counts, not the lasting offset);
               window positive if any (or, as a variant, most) samples are corrupted.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import fftconvolve

from .data import WIN, stage_type
from .decompose import bateman, h_max, t_peak

UP = 8  # fine grid = 8 x 8 Hz = 64 Hz


@dataclass
class Condition:
    """One cell of the grid (docs/METHOD.md section 8)."""

    c: float = 0.5  # share of artifacts caused by movement
    kappa: float = 1.0  # amplitude ~ (I_e / median I)^kappa
    delta_max: float = 1.0  # maximum delay after motion onset, s
    pi_r: float = 0.0  # share of motion events followed by a real response
    classes: tuple = ("C1", "C2", "C3", "C5")  # composition
    class_weights: tuple | None = None  # None = measured UTD shares
    variant: str = "V1"
    prevalence: float = 0.038
    label_rule: str = "any"

    def key(self) -> str:
        cl = "+".join(self.classes)
        return (
            f"c={self.c}|kappa={self.kappa}|dmax={self.delta_max}|pi_r={self.pi_r}|"
            f"cls={cl}|{self.variant}|prev={self.prevalence}|{self.label_rule}"
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


def _segments(stage: np.ndarray):
    starts = np.r_[0, np.where(stage[1:] != stage[:-1])[0] + 1]
    ends = np.r_[starts[1:], stage.size]
    return list(zip(starts, ends, strict=True))


def _block_mean(x: np.ndarray, k: int = UP) -> np.ndarray:
    n = x.size // k
    return x[: n * k].reshape(n, k).mean(1)


def _trapezoid(nf: int, fs_f: float, t0: float, rise: float, plateau: float, fall: float, height: float):
    t = np.arange(nf) / fs_f - t0
    a = np.zeros(nf)
    r = max(rise, 1 / fs_f)
    f = max(fall, 1 / fs_f)
    up = (t >= 0) & (t < r)
    a[up] = t[up] / r
    pl = (t >= r) & (t < r + plateau)
    a[pl] = 1.0
    dn = (t >= r + plateau) & (t < r + plateau + f)
    a[dn] = 1 - (t[dn] - r - plateau) / f
    return height * a


def _step(nf: int, fs_f: float, t0: float, w: float, height: float):
    t = np.arange(nf) / fs_f - t0
    w = max(w, 1 / fs_f)
    a = np.clip(t / w, 0, 1)
    return height * a


class Generator:
    """Builds synthetic participants from estimated parameters."""

    def __init__(self, library: dict, cfg: dict):
        self.lib = library  # {"C1": [params...], ..., "class_shares": {...}, "G": float}
        self.cfg = cfg
        self.g = cfg["generator"]

    # ------------------------------------------------------------ clean part
    def clean(self, motion: dict, phys: dict, rng: np.random.Generator, pi_r: float):
        fs = motion["fs"]
        n = motion["acc"].shape[0]
        nf = n * UP
        fs_f = fs * UP
        stage = motion["stage"]
        tau0, tau1 = phys["tau0"], self.g["tau1"]
        # tonic: one fitted quadruple per protocol stage, stitched without jumps
        s = np.zeros(nf)
        level = None
        for a, b in _segments(stage):
            typ = stage_type(str(stage[a]))
            pool = [t for t in phys["tonic"] if t["stage_type"] == typ] or phys["tonic"]
            q4 = pool[rng.integers(len(pool))]
            tau = np.arange((b - a) * UP) / fs_f
            seg = q4["b"] * tau + q4["A_s"] * np.sin(2 * np.pi * tau / q4["T_t"] + rng.uniform(0, 2 * np.pi))
            start = q4["s0"] if level is None else level
            seg = seg - seg[0] + start
            s[a * UP : b * UP] = seg
            level = seg[-1]
        # responses: spontaneous (Poisson) and motion-evoked (share pi_r)
        imp = np.zeros(nf)
        responses = []
        areas = np.asarray(phys["areas"], float)
        for a, b in _segments(stage):
            typ = stage_type(str(stage[a]))
            lam = _num(phys["rates"].get(typ, {}).get("rate_per_s"))
            if not np.isfinite(lam):  # donor has no clean time of this stage type: use its other rate
                other = [_num(v["rate_per_s"]) for v in phys["rates"].values() if np.isfinite(_num(v["rate_per_s"]))]
                lam = other[0] if other else 0.0
                self.rate_fallbacks = getattr(self, "rate_fallbacks", 0) + 1
            dur = (b - a) / fs
            k = rng.poisson(lam * dur) if np.isfinite(lam) else 0
            for t_k in np.sort(rng.uniform(a / fs, b / fs, k)):
                area = float(areas[rng.integers(areas.size)])
                responses.append({"t": float(t_k), "area": area, "origin": "spontaneous"})
        for ev in motion["events"]:
            if rng.random() < pi_r:
                lat = rng.uniform(*self.g["motion_response_latency_s"])
                t_k = ev["onset"] / fs + lat
                if t_k < n / fs:
                    area = float(areas[rng.integers(areas.size)])
                    responses.append({"t": float(t_k), "area": area, "origin": "motion"})
        hm, tp = h_max(tau0, tau1), t_peak(tau0, tau1)
        for r in responses:
            idx = int(round(r["t"] * fs_f))
            if idx < nf:
                imp[idx] += r["area"]
            r["amp"] = r["area"] * hm
            r["peak_sample"] = int((r["t"] + tp) * fs)
        h = bateman(np.arange(0, self.g["kernel_tau_multiple"] * tau0, 1 / fs_f), tau0, tau1)
        r_f = fftconvolve(imp, h)[:nf]
        clean8 = _block_mean(s + r_f)
        y0 = clean8 + rng.normal(0.0, phys["sigma"], clean8.size)
        return y0, responses

    # ------------------------------------------------------------ artifacts
    def _template(self, cls: str, rng: np.random.Generator):
        pool = self.lib.get(cls) or []
        if pool:
            return dict(pool[rng.integers(len(pool))]), "measured"
        fb = self.g["fallback_templates"].get(cls)
        if fb is None or not self.g["allow_fallback_templates"]:
            raise ValueError(f"no measured {cls} episodes and fallback templates disabled")
        out = {}
        for k, v in fb.items():
            out[k] = rng.uniform(*v) if isinstance(v, list) else v
        return out, "fallback"

    def _artifact_signal(self, cls, prm, t0, scale, n, fs, y0, acc, ev, variant):
        nf, fs_f = n * UP, fs * UP
        if cls == "C1":
            a = _block_mean(_trapezoid(nf, fs_f, t0, prm["tau_r"], 0.0, prm["tau_f"], prm["A"] * scale))
            return a, np.abs(a) >= self.g["label_threshold"]
        if cls in ("C2", "C3"):
            if cls == "C2":
                shape = _trapezoid(nf, fs_f, t0, prm["w"], prm["D"], prm["w"], 1.0)
                size = -prm["depth"] * scale
            else:
                shape = _step(nf, fs_f, t0, prm["w"], 1.0)
                size = prm["delta"] * scale
            shape8 = _block_mean(shape)
            if variant == "V2":
                i0 = min(int(t0 * fs), n - 1)
                lvl = float(np.median(y0[max(0, i0 - int(5 * fs)) : i0 + 1]))
                rho = np.clip(abs(size) / max(lvl, 1e-6), 0.0, self.g["rho_max"])
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
                return piece, moving & (abs(size) >= self.g["label_threshold"])
            return piece, np.abs(piece) >= self.g["label_threshold"]
        if cls == "C4":
            k = int(round(self.cfg["estimation"]["c4"]["smooth_s"] * fs))
            mag = np.sqrt((acc**2).sum(1))
            fast = mag - np.convolve(mag, np.ones(k) / k, "same")
            out = np.zeros(n)
            out[ev["onset"] : ev["end"]] = self.lib["G"] * scale * fast[ev["onset"] : ev["end"]]
            return out, np.abs(out) >= self.g["label_threshold"]
        raise ValueError(cls)

    def generate(self, pid: str, motion: dict, phys: dict, cond: Condition, rng: np.random.Generator) -> SynthRecord:
        fs = motion["fs"]
        acc = motion["acc"]
        n = acc.shape[0]
        y0, responses = self.clean(motion, phys, rng, cond.pi_r)
        thr = self.g["label_threshold"]
        art = np.zeros(n)
        support = np.zeros(n, bool)  # samples labelled as corrupted
        art_by_cls: dict[str, np.ndarray] = {}
        replace = np.zeros(n, bool)
        replace_val = np.zeros(n)
        artifacts = []
        events = motion["events"]
        med_i = float(np.median([e["intensity"] for e in events])) if events else 1.0
        nw = n // WIN
        shares = cond.class_weights or tuple(self.lib["class_shares"].get(c, 0.0) for c in cond.classes)
        shares = np.asarray(shares, float)
        shares = shares / shares.sum() if shares.sum() > 0 else np.full(len(cond.classes), 1 / len(cond.classes))

        def positive_windows():
            m = support | (replace & (np.abs(replace_val - y0) >= thr))
            return m, _window_labels(m, nw, cond.label_rule)

        mask, lab = positive_windows()
        tries = 0
        while lab.mean() < cond.prevalence and tries < self.g["max_artifacts_per_participant"]:
            tries += 1
            cls = cond.classes[rng.choice(len(cond.classes), p=shares)]
            # C4 has no template: its shape is G * (fast part of |a|) inside a motion event
            prm, src = ({}, "model_G") if cls == "C4" else self._template(cls, rng)
            ev = None
            if cls == "C4":
                if not events:
                    break
                ev = events[rng.integers(len(events))]
                t0, caused = ev["onset"] / fs, True
                scale = (ev["intensity"] / med_i) ** cond.kappa
            elif events and rng.random() < cond.c:
                ev = events[rng.integers(len(events))]
                t0 = ev["onset"] / fs + rng.uniform(0, cond.delta_max)
                caused = True
                scale = (ev["intensity"] / med_i) ** cond.kappa
            else:
                t0, caused, scale = rng.uniform(0, n / fs), False, 1.0
            if cls == "C5":
                a = int(t0 * fs)
                b = min(n, a + max(1, int(round(prm["D"] * fs))))
                replace[a:b] = True
                replace_val[a:b] = 0.0 if prm.get("kind", "zero") == "zero" else self.g["saturation_level"]
            else:
                piece, sup = self._artifact_signal(cls, prm, t0, scale, n, fs, y0, acc, ev, cond.variant)
                art += piece
                support |= sup
                art_by_cls[cls] = art_by_cls.get(cls, np.zeros(n, bool)) | sup
            artifacts.append(
                {
                    "cls": cls,
                    "t0": float(t0),
                    "caused_by_motion": caused,
                    "scale": float(scale),
                    "event_onset": None if ev is None else int(ev["onset"]),
                    "template": src,
                    **prm,
                }
            )
            mask, lab = positive_windows()
        y = y0 + art
        y[replace] = replace_val[replace]
        class_windows = {c: _window_labels(a, nw, "any") for c, a in art_by_cls.items()}
        if replace.any():
            class_windows["C5"] = _window_labels(replace & (np.abs(replace_val - y0) >= thr), nw, "any")
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
                "motion_donor": motion["pid"],
                "physiology_donor": phys["pid"],
                "condition": cond.key(),
                "q": q,
                "q_acc": motion.get("q_acc", 0.0),
                "class_windows": class_windows,
            },
        )


def _num(x) -> float:
    """Number from the bank; JSON null (missing estimate) becomes NaN."""
    return float("nan") if x is None else float(x)


def _window_labels(mask: np.ndarray, nw: int, rule: str) -> np.ndarray:
    m = mask[: nw * WIN].reshape(nw, WIN)
    if rule == "any":
        return m.any(1).astype(int)
    if rule == "majority":
        return (m.sum(1) > WIN / 2).astype(int)
    raise ValueError(rule)


def derangement(n: int, rng: np.random.Generator) -> np.ndarray:
    """Random permutation with no fixed point (physiology donor != motion donor)."""
    while True:
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return p
