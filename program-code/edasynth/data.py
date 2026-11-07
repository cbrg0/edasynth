"""This file is the first step of the program after the download: it loads the real
recordings and the expert labels.

It reads the 20 UTD recordings from PhysioNet, with the acceleration and the
electrodermal activity of each person, assigns every moment to a stage of the experiment
and attaches the labels that three experts gave each five-second window. The result is
one object per person, which all later steps use.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .wfdbio import readAnnotations, readRecord

WIN = 40  # A five-second window contains 40 samples at 8 samples per second.


@dataclass
class Recording:
    """One continuous recording of one person."""

    pid: str
    fs: float
    eda: np.ndarray
    acc: np.ndarray  # an array with one row per sample and the columns x, y and z in units of g
    temp: np.ndarray | None = None
    stage: np.ndarray | None = None  # the name of the protocol stage for every sample
    labels3: np.ndarray | None = None  # the labels of the three experts, 0 or 1, with one row per window
    units: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.eda.size)

    @property
    def nWindows(self) -> int:
        return self.n // WIN

    def cropped(self, offset: int) -> "Recording":
        """Return the same recording without its first offset samples; the labels stay
        unchanged.

        This is used when the alignment check found that the expert windows begin offset
        samples after the start of the raw recording.
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
        """Return the label that at least two of the three experts gave each window."""
        if self.labels3 is None:
            raise ValueError("no labels")
        return (self.labels3.sum(axis=1) >= 2).astype(int)


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
    """Map a stage of the protocol to one of the two types, rest or load, used when
    parameters are measured.
    """
    return "rest" if name.startswith("relax") or name == "rest" else "load"


def assignStages(n: int, annotations: list, stages: list[str] = utdStages) -> np.ndarray:
    """Turn the transition annotations of a recording into a stage name for every sample.

    The PhysioNet page says that the annotation file "indicates the time locations and
    labels of the transition states", but it does not explain the codes. The function
    therefore relies on the order of the annotations. When there are as many annotations
    as stages, they are taken as the starts of the stages. When there is one annotation
    fewer, they are taken as the boundaries between stages, and the first stage starts
    at zero. Any other number raises an error, so a mismatch is never accepted silently.
    """
    t = sorted(a.sample for a in annotations)
    if len(t) == len(stages):
        starts = t
    elif len(t) == len(stages) - 1:
        starts = [0, *t]
    else:
        raise ValueError(
            f"{len(t)} annotations; expected {len(stages)} or {len(stages) - 1}. "
            "Check the .atr annotation file and adjust the stage list."
        )
    out = np.empty(n, dtype=object)
    starts = [*starts, n]
    out[: starts[0]] = stages[0]
    for k, name in enumerate(stages):
        out[starts[k] : starts[k + 1]] = name
    return out


def loadUtdLabels(zhangRoot: str | Path) -> pd.DataFrame:
    """Load the expert labels of Zhang et al. for UTD, with one row per five-second window.

    The columns are the group from 0 to 19, the number of the window within the
    participant, the labels of the three experts e1, e2 and e3, and the majority label.
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
    """Load the published feature tables, which are already scaled within each participant:
    24 EDA features and 96 acceleration features.
    """
    root = Path(zhangRoot) / "UTDallas"
    eda = pd.read_csv(root / "UTD_eda.csv", header=None).to_numpy(float)
    acc = pd.read_csv(root / "UTD_acc.csv", header=None).to_numpy(float)
    return eda, acc


def loadUtdParticipant(noneegRoot: str | Path, subject: int, labels: pd.DataFrame | None = None) -> Recording:
    """Load the recording Subject{subject}_AccTempEDA, counted from 1, with its stages and,
    when they are given, its labels.

    The group g in the label files belongs to Subject{g+1}. The module alignment checks
    this correspondence before the data are used.
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


awwFiles = [  # the order used by Zhang et al. in FeatureMatrix_AWW.py: 5 walking days, then 6 resting days
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
    """Load the 11 AWW recordings with their expert labels.

    The header of the raw files lists the columns as Z-axis, Y-axis, X-axis, Battery,
    Celsius and EDA(uS), so the axes are taken by their names and not by their
    positions. The labels are coded as in Zhang et al.: 0, which means not labelled,
    counts as clean, 1 becomes clean (0) and 2 becomes artifact (1).
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
