"""This script extracts the synchronised wrist data of WESAD, as the last section of
study-description/PREREGISTRATION.md describes.

Run it in the folder that holds WESAD.zip. It writes wesad-wrist.npz, which contains for
every participant the wrist acceleration at 32 samples per second in units of 1/64 g,
the wrist EDA at 4 samples per second in microsiemens, and the protocol label, taken at
every 175th sample of the label recorded at 700 samples per second, which gives 4 labels
per second.
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
