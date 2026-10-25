"""edasynth: semi-synthetic test bench for motion-artifact detection in wrist EDA.

The package answers one conditional question: under which assumed links between
wrist movement and electrode contact does a synchronous accelerometer improve the
detection of motion artifacts in electrodermal activity (EDA)?

Real accelerometer recordings (UTD, PhysioNet "noneeg" 1.0.0) are combined with
synthetic EDA whose clean component, artifacts and true labels are known exactly.
See README.md and docs/ for the full method.
"""

__version__ = "1.2.0"
