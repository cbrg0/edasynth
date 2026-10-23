"""Readers: WFDB format 16 with checksum, MIT annotations, Zhang labels, alignment."""

import numpy as np
import pytest
from conftest import NONEEG, ZHANG, needsNoneeg, needsZhang

from edasynth.wfdbio import readAnnotations, readHeader, readRecord


def writeRecord(tmp, digital, gains, baselines, names, unsignedChecksum=False):
    n, k = digital.shape
    digital.astype("<i2").tofile(tmp / "rec.dat")
    lines = [f"rec {k} 8 {n}"]
    for j in range(k):
        col = digital[:, j].astype(np.int64)
        cs = int(col.sum() % 65536)
        if not unsignedChecksum:
            cs = cs - 65536 if cs > 32767 else cs
        lines.append(f"rec.dat 16 {gains[j]}({baselines[j]})/NU 16 0 {col[0]} {cs} 0 {names[j]}")
    (tmp / "rec.hea").write_text("\n".join(lines) + "\n# age: 30\n")


def testFormat16RoundtripAndChecksum(tmp_path):
    rng = np.random.default_rng(0)
    dig = rng.integers(-30000, 30000, size=(100, 2))
    writeRecord(tmp_path, dig, [2000.5, 100.0], [-500, 10], ["EDA", "ax"])
    sig, hdr = readRecord(tmp_path / "rec")
    assert hdr.names == ["EDA", "ax"] and hdr.fs == 8 and hdr.comments == ["age: 30"]
    assert np.allclose(sig[:, 0], (dig[:, 0] + 500) / 2000.5)
    raw = np.fromfile(tmp_path / "rec.dat", dtype="<i2")
    raw[10] += 1
    raw.tofile(tmp_path / "rec.dat")
    with pytest.raises(ValueError, match="checksum"):
        readRecord(tmp_path / "rec")


def testChecksumAcceptsUnsignedHeaderValue(tmp_path):
    # PhysioNet noneeg headers store checksums as unsigned 16-bit numbers (e.g. 57979 for -7557)
    dig = np.full((10, 1), 30000)  # sum 300000 -> 300000 mod 65536 = 37856, signed -27680
    writeRecord(tmp_path, dig, [1.0], [0], ["ay"], unsignedChecksum=True)
    assert " 37856 " in (tmp_path / "rec.hea").read_text()
    sig, hdr = readRecord(tmp_path / "rec")
    assert sig.shape == (10, 1)
    raw = np.fromfile(tmp_path / "rec.dat", dtype="<i2")
    raw[3] += 1
    raw.tofile(tmp_path / "rec.dat")
    with pytest.raises(ValueError, match="checksum"):
        readRecord(tmp_path / "rec")


def testHeaderLineOfUtdParses():
    import tempfile
    from pathlib import Path

    txt = (
        "Subject1_AccTempEDA 5 8 18343\n"
        "Subject1_AccTempEDA.dat 16 19679.8798799(-21943)/NU 16 0 -6986 17497 0 ax\n"
        "Subject1_AccTempEDA.dat 16 24093.3823529(11082)/NU 16 0 -9879 57979 0 ay\n"
        "Subject1_AccTempEDA.dat 16 23239.0070922(11851)/NU 16 0 9527 2729 0 az\n"
        "Subject1_AccTempEDA.dat 16 17245.7894737(-551865)/degC 16 0 -32766 37954 0 temp\n"
        "Subject1_AccTempEDA.dat 16 25361.4551084(-34770)/NU 16 0 -32664 31215 0 EDA\n"
    )
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "s.hea"
        p.write_text(txt)
        h = readHeader(p)
    assert h.nSamples == 18343 and h.names == ["ax", "ay", "az", "temp", "EDA"]
    assert h.signals[4].gain == pytest.approx(25361.4551084) and h.signals[4].baseline == -34770
    assert h.units["EDA"] == "NU" and h.units["temp"] == "degC"


def testMitAnnotationReader(tmp_path):
    words = []

    def ann(code, dt):
        words.append((code << 10) | dt)

    ann(1, 100)
    text = b"Relax"
    words.append((63 << 10) | len(text))
    padded = text + b"\x00" * (len(text) % 2)
    words.extend(np.frombuffer(padded, dtype="<u2").tolist())
    words.append((59 << 10))  # SKIP, then 32-bit increment high word first
    big = 5000
    words.extend([big >> 16, big & 0xFFFF])
    ann(1, 20)
    words.append(0)
    np.asarray(words, dtype="<u2").tofile(tmp_path / "x.atr")
    a = readAnnotations(tmp_path / "x.atr")
    assert [x.sample for x in a] == [100, 5120]
    assert a[0].aux == "Relax"


def testAnnotationReaderMatchesWfdbOnRealData():
    wfdb = pytest.importorskip("wfdb")
    if not (NONEEG / "Subject1_AccTempEDA.atr").exists():
        pytest.skip("noneeg not downloaded")
    ours = readAnnotations(NONEEG / "Subject1_AccTempEDA.atr")
    ref = wfdb.rdann(str(NONEEG / "Subject1_AccTempEDA"), "atr")
    assert [a.sample for a in ours] == list(ref.sample)


@needsNoneeg
def testAllUtdRecordsPassChecksums():
    for s in range(1, 21):
        sig, hdr = readRecord(NONEEG / f"Subject{s}_AccTempEDA")
        assert sig.shape == (hdr.nSamples, 5)


@needsNoneeg
def testRecordReaderMatchesWfdb():
    wfdb = pytest.importorskip("wfdb")
    sig, unused = readRecord(NONEEG / "Subject1_AccTempEDA")
    ref = wfdb.rdrecord(str(NONEEG / "Subject1_AccTempEDA"))
    assert np.allclose(sig, ref.p_signal, atol=1e-9)


@needsZhang
def testUtdLabelsCounts():
    from edasynth.data import loadUtdLabels

    lab = loadUtdLabels(ZHANG)
    assert len(lab) == 9211
    assert int(lab.majority.sum()) == 349
    assert lab.group.nunique() == 20
    assert (lab.groupby("group").size() == 458).iloc[0]


def testAlignmentDetectsAShift():
    from edasynth.alignment import checkParticipant
    from edasynth.features import stats4, windows

    rng = np.random.default_rng(0)
    eda = np.cumsum(rng.standard_normal(40 * 300)) * 0.01 + 5
    stored = stats4(windows(eda))[:280]
    res = checkParticipant(eda, stored, 0.999)
    assert res["accepted"] and res["shiftWindows"] == 0
    res = checkParticipant(np.r_[np.zeros(80), eda], stored, 0.999)
    assert res["accepted"] and res["shiftWindows"] * 40 + res["shiftSamples"] == 80


def testAlignmentOffsetIsAppliedByCropping():
    from edasynth.alignment import checkParticipant
    from edasynth.data import Recording
    from edasynth.features import stats4, windows

    rng = np.random.default_rng(1)
    eda = np.cumsum(rng.standard_normal(40 * 300)) * 0.01 + 5
    stored = stats4(windows(eda))[:280]
    padded = np.r_[np.full(83, 5.0), eda]
    res = checkParticipant(padded, stored, 0.999)
    assert res["accepted"] and res["offsetSamples"] == 83
    r = Recording("S01", 8.0, padded, np.zeros((padded.size, 3)), stage=np.full(padded.size, "relax1", dtype=object))
    c = r.cropped(res["offsetSamples"])
    assert np.array_equal(c.eda[: eda.size], eda)
    assert c.acc.shape[0] == c.eda.size and c.stage.size == c.eda.size
    assert c.meta["alignmentOffsetSamples"] == 83


def testAlignmentRejectsWhenOffsetLeavesTooFewWindows():
    from edasynth.alignment import checkParticipant
    from edasynth.features import stats4, windows

    rng = np.random.default_rng(2)
    eda = np.cumsum(rng.standard_normal(40 * 100)) * 0.01 + 5
    stored = stats4(windows(eda))
    res = checkParticipant(np.r_[np.full(40, 5.0), eda[:-20]], stored, 0.999)
    assert res["nWindowsAfterOffset"] < stored.shape[0]
    assert not res["accepted"]
