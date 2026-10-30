"""Extract the synchronised wrist data of WESAD (docs/EXTERNAL-CHECKS.md, preparation).

Run in the folder that holds WESAD.zip; writes wesad-wrist.npz with, per subject, the
wrist acceleration (32 Hz, 1/64 g), the wrist EDA (4 Hz, µS) and the protocol label
taken every 175th sample of the 700 Hz label (4 Hz).
"""

import pickle
import zipfile

import numpy as np

out = {}
with zipfile.ZipFile("WESAD.zip") as z:
    names = sorted(n for n in z.namelist() if n.endswith(".pkl"))
    for n in names:
        with z.open(n) as f:
            d = pickle.load(f, encoding="latin1")  # noqa: S301 - the published WESAD files
        s = d["subject"]
        w = d["signal"]["wrist"]
        out[f"{s}__acc"] = np.asarray(w["ACC"], np.float32)
        out[f"{s}__eda"] = np.asarray(w["EDA"], np.float32).ravel()
        out[f"{s}__label4hz"] = np.asarray(d["label"])[::175].astype(np.int8)
        print(s, w["ACC"].shape, w["EDA"].shape, len(d["label"]), flush=True)
        del d
np.savez_compressed("wesad-wrist.npz", **out)
print("done", len(names))
