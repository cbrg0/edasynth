"""This file stores the measurements of the previous step in one place, called the bank.

It gathers what estimation.py measured for each participant, groups the artifact shapes
by type and writes everything to two files in result-tables. The generator later reads
only this bank and never the raw recordings.
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np

from .data import Recording
from .estimation import estimateParticipant, quantisationStep, templateParameters

physKeys = (
    "pid",
    "tau0",
    "rates",
    "areas",
    "tonic",
    "sigma",
    "sigmaSource",
    "residualSource",
    "residualSd",
    "shareBelowMinLevel",
    "validDonor",
    "q",
)


def buildLibrary(estimates: list[dict], cfg: dict, log=None) -> dict:
    """Build the artifact templates, the shares of the artifact classes and the walking
    gain from the measurements of all participants.

    Episodes that start below the lowest valid EDA level (rule 1 of Kleckner et al.,
    2018) have no meaningful size relative to the level, so they are left out of the
    templates and of the shares. The realism checks also call this function with half of
    the participants, as section 7 of study-description/METHOD.md describes.
    """
    floor = cfg["estimation"]["episodes"]["levelFloor"]
    allEpisodes = [ep for e in estimates for ep in e["episodes"]]
    episodes = [ep for ep in allEpisodes if ep["level"] >= floor]
    if log and len(episodes) < len(allEpisodes):
        log(f"{len(allEpisodes) - len(episodes)} of {len(allEpisodes)} episodes start below {floor} uS; not used")
    counts = Counter(ep["cls"] for ep in episodes)
    classified = {c: n for c, n in counts.items() if c.startswith("C")}
    tot = sum(v for c, v in classified.items() if c != "C4") or 1
    library = {
        c: [templateParameters(ep, floor) for ep in episodes if ep["cls"] == c] for c in ("C1", "C2", "C3", "C5")
    }
    library["classShares"] = {c: classified.get(c, 0) / tot for c in ("C1", "C2", "C3", "C5")}
    gains = [e["c4"]["G"] for e in estimates if np.isfinite(e["c4"]["G"])]
    library["G"] = float(np.median(gains)) if gains else float("nan")
    library["episodeCounts"] = dict(counts)
    library["episodesBelowFloor"] = len(allEpisodes) - len(episodes)
    return library


def gitCommit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001 - recording the program version is helpful but not essential
        return "unknown"


def buildBank(recordings: list[Recording], cfg: dict, source: str, out: str | Path, log=print, nJobs: int = 1) -> dict:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    log(f"estimating {len(recordings)} recordings on {nJobs} job(s) ...")
    if nJobs == 1:
        estimates = []
        for rec in recordings:
            log(f"estimating {rec.pid} ...")
            estimates.append(estimateParticipant(rec, cfg))
    else:
        from joblib import Parallel, delayed

        estimates = Parallel(n_jobs=nJobs)(delayed(estimateParticipant)(rec, cfg) for rec in recordings)
    arrays = {}
    for e in estimates:  # The noise blocks are stored in the .npz file and not in the JSON file.
        arrays[f"resid{e['pid']}"] = np.asarray(e.pop("residual"), float)
    library = buildLibrary(estimates, cfg, log)
    phys = [{k: e[k] for k in physKeys} for e in estimates]
    excluded = [p["pid"] for p in phys if not p["validDonor"]]
    if excluded:
        log(f"not used as physiology donors (EDA below the minimum level): {', '.join(excluded)}")
    motion = []
    for rec, e in zip(recordings, estimates, strict=True):
        motion.append(
            {
                "pid": rec.pid,
                "fs": rec.fs,
                "events": e["motionEvents"],
                "qAcc": quantisationStep(rec.acc[:, 0]),
                "theta": e["motionTheta"],
            }
        )
        arrays[f"acc{rec.pid}"] = rec.acc
        arrays[f"stage{rec.pid}"] = rec.stage.astype(str)
    meta = {
        "source": source,
        "gitCommit": gitCommit(),
        "config": cfg.get("configPath"),
        "phys": phys,
        "motion": motion,
        "library": library,
        "estimates": estimates,
    }
    out.with_suffix(".json").write_text(toJson(meta))
    np.savez_compressed(out.with_suffix(".npz"), **arrays)
    return meta


def toJson(obj) -> str:
    """Write JSON in which missing or infinite numbers become null, because standard JSON
    has no NaN.
    """
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

    return json.dumps(clean(obj), indent=1, default=jsonDefault, allow_nan=False)


def jsonDefault(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))
