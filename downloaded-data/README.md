# Data

The data are not stored in the repository. The command `make data` downloads them from
their sources into this folder and compares every file with a list of checksums. A
checksum is a short code computed from the content of a file, so a matching checksum
shows that everyone works with exactly the same files.

| Data | Where they come from | Licence |
|---|---|---|
| Wristband recordings of 20 people, called UTD in this project. | PhysioNet, the Non-EEG Dataset for Assessment of Neurological Status, version 1.0.0, by Birjandtalab et al. (2016, 2017). | Open Data Commons Attribution License 1.0 |
| The artifact labels that three experts gave these recordings. | The repository of Zhang et al. (2017) at github.com/IdeasLabUT/EDA-Artifact-Detection, at commit `7c4bd46`. | BSD licence |

PhysioNet does not state the unit of the EDA values. They are treated as microsiemens,
because other recordings from the same device, the Affectiva Q Sensor, which Zhang et
al. (2017) published, give EDA in microsiemens.

The data sets EDABE and WESAD are used only for the checks on independent data. They
are downloaded by hand from the pages of their authors and prepared with the scripts
in the folder `data-download-and-external-checks/`.
