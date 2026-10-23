"""Turning per-participant estimates into the generator's inputs ("bank").

Output: ``<name>.json`` (all estimated numbers with provenance) and ``<name>.npz``
(accelerometer arrays and stage labels of the motion donors).
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np

from .data import Recording
from .estimation import estimate_participant, quantisation_step, template_parameters


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance is best effort
        return "unknown"


def build_bank(
    recordings: list[Recording], cfg: dict, source: str, out: str | Path, log=print, n_jobs: int = 1
) -> dict:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    log(f"estimating {len(recordings)} recordings on {n_jobs} job(s) ...")
    if n_jobs == 1:
        estimates = []
        for rec in recordings:
            log(f"estimating {rec.pid} ...")
            estimates.append(estimate_participant(rec, cfg))
    else:
        from joblib import Parallel, delayed

        estimates = Parallel(n_jobs=n_jobs)(delayed(estimate_participant)(rec, cfg) for rec in recordings)
    episodes = [ep for e in estimates for ep in e["episodes"]]
    counts = Counter(ep["cls"] for ep in episodes)
    classified = {c: n for c, n in counts.items() if c.startswith("C")}
    tot = sum(v for c, v in classified.items() if c != "C4") or 1
    library = {c: [template_parameters(ep) for ep in episodes if ep["cls"] == c] for c in ("C1", "C2", "C3", "C5")}
    library["class_shares"] = {c: classified.get(c, 0) / tot for c in ("C1", "C2", "C3", "C5")}
    gains = [e["c4"]["G"] for e in estimates if np.isfinite(e["c4"]["G"])]
    library["G"] = float(np.median(gains)) if gains else float("nan")
    library["episode_counts"] = dict(counts)
    phys = [
        {k: e[k] for k in ("pid", "tau0", "rates", "areas", "tonic", "sigma", "sigma_source", "q")} for e in estimates
    ]
    motion, arrays = [], {}
    for rec, e in zip(recordings, estimates, strict=True):
        motion.append(
            {
                "pid": rec.pid,
                "fs": rec.fs,
                "events": e["motion_events"],
                "q_acc": quantisation_step(rec.acc[:, 0]),
                "theta": e["motion_theta"],
            }
        )
        arrays[f"acc_{rec.pid}"] = rec.acc
        arrays[f"stage_{rec.pid}"] = rec.stage.astype(str)
    meta = {
        "source": source,
        "git_commit": git_commit(),
        "config": cfg.get("_config_path"),
        "phys": phys,
        "motion": motion,
        "library": library,
        "estimates": estimates,
    }
    out.with_suffix(".json").write_text(to_json(meta))
    np.savez_compressed(out.with_suffix(".npz"), **arrays)
    return meta


def to_json(obj) -> str:
    """JSON with NaN/inf written as null (standard JSON has no NaN)."""
    import math

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, np.ndarray):
            return clean(o.tolist())
        if isinstance(o, (float, np.floating)) and not math.isfinite(float(o)):
            return None
        return o

    return json.dumps(clean(obj), indent=1, default=_default, allow_nan=False)


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))
