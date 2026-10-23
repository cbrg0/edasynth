# Data

This folder receives the two public data sets used by the study. They are not stored
in the repository; `make data` (script `scripts/download-data.sh`) downloads them at
fixed versions and checks every file against a list of SHA-256 checksums. Only that
checksum list is part of the repository.

Files under `raw/` and `external/` are never edited by hand or by code. Every
processing step writes its output to `results/`.

## Contents after `make data`

| Folder | Source | Version | Licence | Used for |
|---|---|---|---|---|
| `raw/noneeg-1.0.0/` | PhysioNet, *Non-EEG Dataset for Assessment of Neurological Status* (Birjandtalab et al., 2016, 2017) | 1.0.0, DOI 10.13026/C26Q2Z | Open Data Commons Attribution License v1.0 | `SubjectN_AccTempEDA.{hea,dat,atr}` for N = 1…20: three acceleration axes, skin temperature and EDA at 8 Hz, and the protocol-stage annotations |
| `external/eda-artifact-detection/` | Zhang et al. (2017), https://github.com/IdeasLabUT/EDA-Artifact-Detection | commit `7c4bd467bec9099669ddf3581b00a9a843a2f136` | BSD | expert artifact labels for UTD (`UTDallas/*.csv`), the published feature matrices, and the raw AWW recordings with their labels |
| `PROVENANCE.txt` | written by the download script | — | — | download date and sources |

The PhysioNet archive also contains `SubjectN_SpO2HR` records (heart rate and blood
oxygen). They are downloaded because they are part of the archive but are not used.

## Integrity checks

- PhysioNet files are checked against the `SHA256SUMS.txt` file that ships with the
  data.
- The Zhang et al. files are checked against
  `checksums/zhang-eda-artifact-detection-7c4bd46.sha256` (43 files). This list was
  computed from a copy of the repository at the pinned commit.

`python -m edasynth check-data` repeats both checks at any time.

## Facts about the files

These facts were read from the files and are used by the code:

- The header of `Subject1_AccTempEDA` lists 5 signals at 8 Hz, 18 343 samples, stored
  in WFDB format 16. The units are given as `NU` (no unit) for `ax`, `ay`, `az` and
  `EDA`, and `degC` for `temp`. PhysioNet therefore does not state the unit of EDA.
  The values are treated as µS, which is consistent with the device (Affectiva Q
  Sensor, according to Zhang et al., 2017) and with the `EDA(uS)` column of the AWW
  files from the same device, but this is an assumption until `make inspect` has been
  run.
- The label files contain 9 211 windows of 5 s in 20 groups numbered 0–19; 349
  windows (3.79 %) are artifacts by majority of three experts. Group 0 has
  458 windows, which equals the number of complete 40-sample windows in
  `Subject1` (18 343 // 40).
- PhysioNet describes the annotation file as marking the time locations and labels
  of the transition states, without documenting the codes. The stage parser therefore
  relies on the order of the annotations and stops with an error if their number is
  unexpected.
