#!/usr/bin/env python3

import json
import platform

import numpy as np
import pandas as pd
import sklearn
import yaml
from sklearn.cluster import DBSCAN, KMeans


def parse_eigenvec(path):
    """Parse a PLINK2 .eigenvec file into (sample_ids: pd.Series, pcs: np.ndarray).

    Accepts the FID/IID and IID-only header layouts PLINK2 emits, plus the
    leading '#' on the header line. Sample IDs are read from the IID column.
    """
    df = pd.read_csv(path, sep=r"\\s+", engine="python")
    df.columns = [c.lstrip("#") for c in df.columns]
    cols_upper = [c.upper() for c in df.columns]
    if cols_upper[:2] == ["FID", "IID"]:
        id_cols = df.columns[:2]
    elif cols_upper[:1] == ["IID"]:
        id_cols = df.columns[:1]
    else:
        raise ValueError(f"eigenvec file missing IID header: {list(df.columns)}")
    sample_ids = df["IID"].astype(str)
    pcs = df.drop(columns=id_cols).to_numpy(dtype=float)
    return sample_ids, pcs


def main():
    eigenvec = "$eigenvec"
    algorithm = "$algorithm"
    n_clusters = int("$n_clusters")
    dbscan_eps = float("$dbscan_eps")
    dbscan_min_samples = int("$dbscan_min_samples")
    prefix = "${task.ext.prefix ?: meta.id}"

    sample_ids, x = parse_eigenvec(eigenvec)

    if algorithm == "kmeans":
        model = KMeans(n_clusters=n_clusters, init="random", n_init=100, random_state=42)
        labels = model.fit_predict(x)
        info = {
            "algorithm": "kmeans",
            "k": n_clusters,
            "inertia": float(model.inertia_),
        }
    elif algorithm == "dbscan":
        model = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples)
        labels = model.fit_predict(x)
        info = {
            "algorithm": "dbscan",
            "eps": dbscan_eps,
            "min_samples": dbscan_min_samples,
            "n_clusters_found": len(set(labels) - {-1}),
            "n_noise": int(np.sum(labels == -1)),
        }
    else:
        raise ValueError(f"Unknown algorithm '{algorithm}' (expected 'kmeans' or 'dbscan')")

    info |= {"n_samples": int(x.shape[0]), "n_features": int(x.shape[1])}

    pd.DataFrame({"sample_id": sample_ids, "cluster": labels}).to_csv(f"{prefix}.clusters.csv", index=False)
    with open(f"{prefix}.clustering_info.json", "w") as fh:
        json.dump(info, fh, indent=2)

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "scikit-learn": sklearn.__version__,
        }
    }
    with open("versions.yml", "w") as fh:
        yaml.dump(versions, fh, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    main()
