"""Minimal reader for the PhysioNet WFDB files of the UTD ("noneeg") dataset.

Only what the dataset needs is implemented: signal format 16 (little-endian
16-bit integers, channels interleaved) and MIT-format annotation files (.atr).
The reader verifies every signal against the checksum and initial value stored in
the header, so a wrong byte order, channel order or gain is detected immediately.
If the official ``wfdb`` package is installed, ``tests/test_io_wfdb.py`` also
compares both readers on the real files.

Physical value = (digital value - baseline) / gain   (WFDB convention).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_SIG_RE = re.compile(
    r"^(?P<file>\S+)\s+(?P<fmt>\d+)\s+"
    r"(?P<gain>[-+0-9.eE]+)(?:\((?P<baseline>-?\d+)\))?(?:/(?P<units>\S+))?\s+"
    r"(?P<res>\d+)\s+(?P<zero>-?\d+)\s+(?P<init>-?\d+)\s+(?P<checksum>-?\d+)\s+"
    r"(?P<block>\d+)\s*(?P<name>.*)$"
)


@dataclass
class SignalSpec:
    """One signal line of a WFDB header."""

    file: str
    fmt: int
    gain: float
    baseline: int
    units: str
    adc_res: int
    adc_zero: int
    init: int
    checksum: int
    name: str


@dataclass
class Header:
    """Parsed WFDB header."""

    record: str
    n_sig: int
    fs: float
    n_samples: int
    signals: list[SignalSpec]
    comments: list[str] = field(default_factory=list)

    @property
    def names(self) -> list[str]:
        return [s.name for s in self.signals]

    @property
    def units(self) -> dict[str, str]:
        return {s.name: s.units for s in self.signals}


def read_header(path: str | Path) -> Header:
    """Parse a ``.hea`` file."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    comments = [ln[1:].strip() for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if ln.strip() and not ln.startswith("#")]
    rec = body[0].split()
    record, n_sig, fs, n_samples = rec[0], int(rec[1]), float(rec[2]), int(rec[3])
    sigs = []
    for ln in body[1 : 1 + n_sig]:
        m = _SIG_RE.match(ln.strip())
        if m is None:
            raise ValueError(f"cannot parse signal line: {ln!r}")
        gain = float(m["gain"])
        baseline = int(m["baseline"]) if m["baseline"] is not None else int(m["zero"])
        sigs.append(
            SignalSpec(
                file=m["file"],
                fmt=int(m["fmt"]),
                gain=gain if gain != 0 else 200.0,  # WFDB default gain
                baseline=baseline,
                units=m["units"] or "mV",
                adc_res=int(m["res"]),
                adc_zero=int(m["zero"]),
                init=int(m["init"]),
                checksum=int(m["checksum"]),
                name=m["name"].strip(),
            )
        )
    return Header(record, n_sig, fs, n_samples, sigs, comments)


def _checksum16(digital: np.ndarray) -> int:
    cs = int(np.sum(digital.astype(np.int64)) % 65536)
    return cs - 65536 if cs > 32767 else cs


def read_record(record_path: str | Path, verify: bool = True) -> tuple[np.ndarray, Header]:
    """Read a format-16 record.

    Args:
        record_path: path without extension, e.g. ``data/raw/noneeg/Subject1_AccTempEDA``.
        verify: check initial values and checksums against the header.

    Returns:
        ``(signals, header)`` where ``signals`` has shape (n_samples, n_sig) in physical units.
    """
    record_path = Path(record_path)
    hdr = read_header(record_path.with_suffix(".hea"))
    files = {s.file for s in hdr.signals}
    if len(files) != 1 or any(s.fmt != 16 for s in hdr.signals):
        raise NotImplementedError("only single-file format-16 records are supported")
    dat = record_path.parent / files.pop()
    raw = np.fromfile(dat, dtype="<i2")
    if raw.size != hdr.n_samples * hdr.n_sig:
        raise ValueError(f"{dat.name}: {raw.size} values, header expects {hdr.n_samples * hdr.n_sig}")
    digital = raw.reshape(hdr.n_samples, hdr.n_sig)
    if verify:
        for k, s in enumerate(hdr.signals):
            if digital[0, k] != s.init:
                raise ValueError(f"{s.name}: initial value {digital[0, k]} != header {s.init}")
            cs = _checksum16(digital[:, k])
            if cs != s.checksum:
                raise ValueError(f"{s.name}: checksum {cs} != header {s.checksum}")
    phys = np.empty(digital.shape, dtype=float)
    for k, s in enumerate(hdr.signals):
        phys[:, k] = (digital[:, k].astype(float) - s.baseline) / s.gain
    return phys, hdr


@dataclass
class Annotation:
    """One MIT-format annotation."""

    sample: int
    code: int
    aux: str = ""
    sub: int = 0
    chan: int = 0
    num: int = 0


def read_annotations(path: str | Path) -> list[Annotation]:
    """Read an MIT-format annotation file (e.g. ``.atr``).

    Format: a stream of 16-bit little-endian words; the top 6 bits hold the code
    A, the low 10 bits hold I. A=59 (SKIP) is followed by a 32-bit time increment
    stored as two words, high word first; A=60/61/62 set NUM/SUB/CHAN; A=63 (AUX)
    is followed by I bytes of text (padded to an even length); A=0 with I=0 ends
    the file. Any other A is an annotation type with time increment I.
    """
    buf = np.fromfile(Path(path), dtype="<u2")
    out: list[Annotation] = []
    t = 0
    num = sub = chan = 0
    i = 0
    while i < buf.size:
        w = int(buf[i])
        a, val = w >> 10, w & 0x3FF
        i += 1
        if a == 0 and val == 0:
            break
        if a == 59:  # SKIP
            hi, lo = int(buf[i]), int(buf[i + 1])
            i += 2
            inc = (hi << 16) | lo
            if inc >= 2**31:
                inc -= 2**32
            t += inc
        elif a == 60:
            num = val
        elif a == 61:
            sub = val
        elif a == 62:
            chan = val
        elif a == 63:  # AUX string attached to the previous annotation
            nbytes = val
            nwords = (nbytes + 1) // 2
            raw = buf[i : i + nwords].astype("<u2").tobytes()[:nbytes]
            i += nwords
            if out:
                out[-1].aux = raw.decode("latin-1").rstrip("\x00")
        else:
            t += val
            out.append(Annotation(sample=t, code=a, sub=sub, chan=chan, num=num))
    return out
