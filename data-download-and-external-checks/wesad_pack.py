"""This script prepares the WESAD data set for the external validation.

It takes from the WESAD archive only what the study needs: the wrist acceleration, the
wrist electrodermal activity and the stage of the experiment for every participant. Run
it in the folder with WESAD.zip; it writes wesad-wrist.npz.
"""

import pickle
import zipfile

import numpy as np

out = {}
with zipfile.ZipFile("WESAD.zip") as z:
    names = sorted(n for n in z.namelist() if n.endswith(".pkl"))
    for n in names:
        with z.open(n) as f:
            d = pickle.load(f, encoding="latin1")  # noqa: S301 - these are the published WESAD files
        s = d["subject"]
        w = d["signal"]["wrist"]
        out[f"{s}__acc"] = np.asarray(w["ACC"], np.float32)
        out[f"{s}__eda"] = np.asarray(w["EDA"], np.float32).ravel()
        out[f"{s}__label4hz"] = np.asarray(d["label"])[::175].astype(np.int8)
        print(s, w["ACC"].shape, w["EDA"].shape, len(d["label"]), flush=True)
        del d
np.savez_compressed("wesad-wrist.npz", **out)
print("done", len(names))
