# Data

The data are not stored in the repository. The command `make data` downloads them from
their sources into this folder and checks every file against a list of checksums, so
that everyone works with exactly the same files.

| Data | Source | Licence |
|---|---|---|
| wristband recordings of 20 people (UTD) | PhysioNet, Non-EEG Dataset, version 1.0.0 (Birjandtalab et al., 2016, 2017) | ODC-By 1.0 |
| expert artifact labels for these recordings | Zhang et al. (2017), github.com/IdeasLabUT/EDA-Artifact-Detection, commit `7c4bd46` | BSD |

PhysioNet does not state the unit of EDA; the values are treated as microsiemens, as
on the device used (Affectiva Q Sensor). The data sets EDABE and WESAD, used only for
the independent checks, are downloaded by hand from their authors' pages and prepared
with the scripts in `data-download-and-external-checks/`.
