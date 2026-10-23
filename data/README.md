# Data

Nothing in this folder is committed except the checksum list. `make data`
(`scripts/download_data.sh`) downloads both sources at fixed versions and verifies
every file. Files under `raw/` and `external/` are never edited; every processing
step writes its output to `results/`.

| Folder | Source | Version | Licence | What we use |
|---|---|---|---|---|
| `raw/noneeg-1.0.0/` | PhysioNet, *Non-EEG Dataset for Assessment of Neurological Status* (Birjandtalab et al., 2016, 2017) | 1.0.0, DOI 10.13026/C26Q2Z | Open Data Commons Attribution 1.0 | `SubjectN_AccTempEDA.{hea,dat,atr}` for N = 1..20: 3-axis acceleration, skin temperature and EDA at 8 Hz; protocol-stage annotations |
| `external/eda-artifact-detection/` | Zhang et al. (2017), github.com/IdeasLabUT/EDA-Artifact-Detection | commit `7c4bd467bec9099669ddf3581b00a9a843a2f136` | BSD | expert artifact labels for UTD (`UTDallas/*.csv`), published feature matrices, raw AWW recordings |

Integrity: PhysioNet files are checked against the official `SHA256SUMS.txt`
shipped with the data; the Zhang files against `checksums/zhang-eda-artifact-detection-7c4bd46.sha256`
(computed from the pinned commit).

## Facts verified so far

* Header of `Subject1_AccTempEDA`: 5 signals, 8 Hz, 18 343 samples, format 16.
  Units: `ax, ay, az, EDA` = **NU (no unit given)**, `temp` = degC. The EDA unit is
  therefore not stated by PhysioNet. The representable EDA range of the file
  (about 0.08-2.66 for subject 1) is compatible with microsiemens, and the same
  device writes `EDA(uS)` in the AWW files, but this remains an assumption until
  `make inspect` is run (see `docs/STATUS.md`).
* Labels: 9 211 windows of 5 s, 349 artifacts by majority of three experts
  (3.79 %), 20 groups numbered 0-19; group 0 has 458 windows = floor(18 343 / 40).
* PhysioNet describes the annotation file as marking "the time locations and
  labels of the transition states" without documenting the codes; the stage
  parser therefore relies on their order and refuses any unexpected count.

## Not used

`SubjectN_SpO2HR` records (heart rate, SpO2) are not needed for the question and
are only downloaded because they are part of the archive.
