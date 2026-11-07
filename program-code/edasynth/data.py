"""This file belongs to the first step of the study, preparing the data. It loads the real
recordings and the expert labels so that all later steps work with the same objects.

The recordings are the UTD data set published on PhysioNet (Birjandtalab et al., 2016,
2017): 20 people wore a wristband that recorded wrist acceleration and electrodermal
activity 8 times per second. The labels come from Zhang et al. (2017), in whose study
three experts marked every five-second window as clean or as an artifact. The
five-second window is taken from these authors and from Taylor et al. (2015), who
labelled artifacts in the same way, so the program uses the unit in which the labels
exist. A window counts as an artifact when at least two of the three experts marked it.
This majority rule is a decision of this study; it keeps a window that only one expert
doubted out of the artifacts.

The PhysioNet files mark the transitions between the stages of the experiment but do not
explain their codes. The file therefore reads the stages by their order, and stops with
an error if their number does not fit, so that a wrong reading cannot pass unnoticed.
The outcome is one object per person with the signals, the stages and the labels.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .wfdbio import readAnnotations, readRecord

WIN = 40


@dataclass
class Recording:

    pid: str
    fs: float
    eda: np.ndarray
    acc: np.ndarray
    temp: np.ndarray | None = None
    stage: np.ndarray | None = None
    labels3: np.ndarray | None = None
    units: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.eda.size)

    @property
    def nWindows(self) -> int:
        return self.n // WIN

    def cropped(self, offset: int) -> "Recording":
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
    return "rest" if name.startswith("relax") or name == "rest" else "load"


def assignStages(n: int, annotations: list, stages: list[str] = utdStages) -> np.ndarray:
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
    root = Path(zhangRoot) / "UTDallas"
    eda = pd.read_csv(root / "UTD_eda.csv", header=None).to_numpy(float)
    acc = pd.read_csv(root / "UTD_acc.csv", header=None).to_numpy(float)
    return eda, acc


def loadUtdParticipant(noneegRoot: str | Path, subject: int, labels: pd.DataFrame | None = None) -> Recording:
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


awwFiles = [
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
