#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp"
os.environ["NUMBA_CACHE_DIR"] = "./tmp"

import anndata as ad
import doubletdetection
import yaml

adata = ad.read_h5ad("${h5ad}")

clf = doubletdetection.BoostClassifier(n_jobs=int("${task.cpus}"))
doublets = clf.fit(adata.X).predict()
scores = clf.doublet_score()

adata.obs["doublet"] = [label == 1 for label in doublets]
adata.obs["doublet_score"] = scores

adata.write_h5ad("${prefix}.h5ad")

df = adata.obs[["doublet"]]
df.columns = ["${prefix}"]
df.to_pickle("${prefix}.pkl")

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "doubletdetection": doubletdetection.__version__,
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
