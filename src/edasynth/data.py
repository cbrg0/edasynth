"""Loading the two real data sources.

UTD (PhysioNet ``noneeg`` 1.0.0, Birjandtalab et al., 2016/2017)
    20 participants, wrist Affectiva Q Sensor, 8 Hz, channels ax, ay, az, temp, EDA.
    Raw signals come from PhysioNet; expert artifact labels come from the
    repository of Zhang et al. (2017), one label per non-overlapping 5-s window.

AWW (Alan Walks Wales, distributed with Zhang et al., 2017)
    One person, everyday recordings with the same device. Used only for checking
    the feature code against the published feature matrices and as a development
    source when UTD raw data are not available (``smoke`` runs).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .wfdbio import readAnnotations, readRecord

WIN = 40  # samples per 5-s window at 8 Hz


@dataclass
class Recording:
    """One continuous recording of one person."""

    pid: str
    fs: float
    eda: np.ndarray
    acc: np.ndarray  # (n, 3) columns x, y, z in g
    temp: np.ndarray | None = None
    stage: np.ndarray | None = None  # per-sample stage name
    labels3: np.ndarray | None = None  # (nWindows, 3) expert labels 0/1
    units: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.eda.size)

    @property
    def nWindows(self) -> int:
        return self.n // WIN

    def cropped(self, offset: int) -> "Recording":
        """The same recording without its first ``offset`` samples (labels unchanged).

        Used when the alignment check found that the expert windows start
        ``offset`` samples after the beginning of the raw record.
        """
        if offset <= 0:
            return self
        cut = lambda a: None if a is None else a[offset:]  # noqa: E731
        return Recording(
            pid=self.pid,
            fs=self.fs,
            eda=self.eda[offset:],
            acc=self.acc[offset:],
            temp=cut(self.temp),
            stage=cut(self.stage),
            labels3=self.labels3,
            units=self.units,
            meta={**self.meta, "alignmentOffsetSamples": int(offset)},
        )

    def majority(self) -> np.ndarray:
        """Majority label (at least two of three experts) per window."""
        if self.labels3 is None:
            raise ValueError("no labels")
        return (self.labels3.sum(axis=1) >= 2).astype(int)


# --------------------------------------------------------------------------- UTD
utdStages = [
    "relax1",
    "physical",
    "relax2",
    "instructions",
    "cognitive",
    "relax3",
    "emotional",
    "relax4",
]


def stageType(name: str) -> str:
    """Map a protocol stage to the two types used for parameter estimation."""
    return "rest" if name.startswith("relax") or name == "rest" else "load"


def assignStages(n: int, annotations: list, stages: list[str] = utdStages) -> np.ndarray:
    """Turn transition annotations into a per-sample stage label.

    The PhysioNet page states that the annotation file "indicates the time
    locations and labels of the transition states" but does not document the
    codes. We therefore rely on the order of annotations: with len(stages)
    annotations they are taken as stage starts; with len(stages) - 1 they are
    taken as boundaries (first stage starts at 0). Any other count raises, so a
    mismatch is never silently accepted.
    """
    t = sorted(a.sample for a in annotations)
    if len(t) == len(stages):
        starts = t
    elif len(t) == len(stages) - 1:
        starts = [0, *t]
    else:
        raise ValueError(
            f"{len(t)} annotations; expected {len(stages)} or {len(stages) - 1}. "
            "Inspect the .atr file (python -m edasynth inspect) and adjust the stage list."
        )
    out = np.empty(n, dtype=object)
    starts = [*starts, n]
    out[: starts[0]] = stages[0]
    for k, name in enumerate(stages):
        out[starts[k] : starts[k + 1]] = name
    return out


def loadUtdLabels(zhangRoot: str | Path) -> pd.DataFrame:
    """Expert labels of Zhang et al. for UTD: one row per 5-s window.

    Columns: group (0-19), window (index within participant), e1, e2, e3, majority.
    """
    root = Path(zhangRoot) / "UTDallas"
    groups = pd.read_csv(root / "UTD_groups.csv", header=None)[0].astype(int).to_numpy()
    all3 = pd.read_csv(root / "UTD_label_All3.csv", header=None).to_numpy().astype(int)
    maj = pd.read_csv(root / "UTD_label.csv", header=None)[0].astype(int).to_numpy()
    df = pd.DataFrame({"group": groups, "e1": all3[:, 0], "e2": all3[:, 1], "e3": all3[:, 2]})
    df["majority"] = maj
    df["window"] = df.groupby("group").cumcount()
    if not np.array_equal((all3.sum(1) >= 2).astype(int), maj):
        raise ValueError("UTD_label.csv is not the majority of UTD_label_All3.csv")
    return df


def loadUtdStoredFeatures(zhangRoot: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Published feature matrices (already z-scored within participant): EDA 24, ACC 96."""
    root = Path(zhangRoot) / "UTDallas"
    eda = pd.read_csv(root / "UTD_eda.csv", header=None).to_numpy(float)
    acc = pd.read_csv(root / "UTD_acc.csv", header=None).to_numpy(float)
    return eda, acc


def loadUtdParticipant(noneegRoot: str | Path, subject: int, labels: pd.DataFrame | None = None) -> Recording:
    """Load ``Subject{subject}_AccTempEDA`` (1-based) with stages and, if given, labels.

    Group ``g`` in the label files corresponds to ``Subject{g+1}``; this mapping
    is verified by :mod:`edasynth.alignment` before any use.
    """
    rec = Path(noneegRoot) / f"Subject{subject}_AccTempEDA"
    sig, hdr = readRecord(rec)
    idx = {name: k for k, name in enumerate(hdr.names)}
    acc = sig[:, [idx["ax"], idx["ay"], idx["az"]]]
    r = Recording(
        pid=f"S{subject:02d}",
        fs=hdr.fs,
        eda=sig[:, idx["EDA"]],
        acc=acc,
        temp=sig[:, idx["temp"]],
        units=hdr.units,
        meta={"headerComments": hdr.comments, "nSamples": hdr.nSamples},
    )
    try:
        ann = readAnnotations(rec.with_suffix(".atr"))
        r.meta["annotations"] = [(a.sample, a.code, a.aux) for a in ann]
        r.stage = assignStages(r.n, ann)
    except (FileNotFoundError, ValueError) as err:
        r.meta["stageError"] = str(err)
    if labels is not None:
        g = labels[labels.group == subject - 1]
        r.labels3 = g[["e1", "e2", "e3"]].to_numpy()
    return r


# --------------------------------------------------------------------------- AWW
awwFiles = [  # order used by Zhang et al. (FeatureMatrix_AWW.py): 5 walking, 6 resting
    "2013_04_24_1hour_walk&chat&drink",
    "2013_05_06 - 1hour_walking",
    "2013_05_28_1hour_walking",
    "2013_06_04_1hour_walking",
    "2013_06_10_1hour_walking",
    "2013_04_24_1hour_chat&drink_Pub",
    "2013_05_14_40mins_eating",
    "2013_05_19_40mins_eat&drink",
    "2013_06_03_1hour_lunch",
    "2013_06_11_40mins_eat&readnewspaper",
    "2013_07_10_1hour_Dinner",
]


def awwColumns(path: Path) -> list[str]:
    with path.open(encoding="latin-1") as fh:
        rows = [next(csv.reader(fh)) for iteration in range(8)]
    header = rows[6][0]
    return [c.strip() for c in header.split("|")]


def loadAww(zhangRoot: str | Path) -> list[Recording]:
    """Load the 11 AWW recordings with expert labels.

    The raw file header lists columns as ``Z-axis | Y-axis | X-axis | Battery |
    Celsius | EDA(uS)``; axes are therefore taken by name, not by position.
    Label coding follows Zhang et al.: 0 (not labelled) is treated as clean,
    then 1 -> clean (0) and 2 -> artifact (1).
    """
    root = Path(zhangRoot) / "AlanWalksWales" / "Raw"
    out = []
    for name in awwFiles:
        f = root / f"{name}.csv"
        cols = awwColumns(f)
        data = np.loadtxt(f, delimiter=",", skiprows=8)
        col = {c.split("-")[0].split("(")[0].strip().lower(): k for k, c in enumerate(cols)}
        lab = np.loadtxt(root / f"{name}_Epochs.csv", delimiter=",", skiprows=1, usecols=(3, 4, 5))
        lab[lab == 0] = 1
        lab = (lab - 1).astype(int)
        n = len(lab) * WIN
        acc = data[:n, [col["x"], col["y"], col["z"]]]
        out.append(
            Recording(
                pid=f"AWW{len(out):02d}",
                fs=8.0,
                eda=data[:n, col["eda"]],
                acc=acc,
                temp=data[:n, [k for c, k in col.items() if "celsius" in c][0]],
                stage=np.full(n, "walk" if "walk" in name else "rest", dtype=object),
                labels3=lab,
                units={"EDA": "uS", "acc": "g"},
                meta={"file": name, "columns": cols},
            )
        )
    return out
