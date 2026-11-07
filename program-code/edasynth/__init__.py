"""The edasynth program studies whether a wrist accelerometer helps to detect motion
artifacts in electrodermal activity (EDA).

The answer depends on how movement and the contact of the electrodes are linked, and
this link is unknown. The program therefore combines real accelerometer recordings of
the UTD data set (PhysioNet, "noneeg", version 1.0.0) with synthetic EDA signals whose
clean part, artifacts and true labels are known exactly, and measures the benefit of the
accelerometer for many assumed links. The method is described in README.md and in the
folder study-description.
"""

__version__ = "1.4.1"
