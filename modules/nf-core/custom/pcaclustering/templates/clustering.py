#!/usr/bin/env python3

# Copyright (c) nf-core
# This software is licensed under the MIT License.
# SPDX-License-Identifier: MIT

import json
import platform

import numpy as np
import pandas as pd
import sklearn
import yaml
from sklearn.cluster import DBSCAN, KMeans


def load_features(path):
    """Read a TSV where the first column is sample IDs and the remaining
    columns are numeric features.

    Returns (sample_ids: pd.Series, features: np.ndarray).
    """
    df = pd.read_csv(path, sep="\\t")
    if df.shape[1] < 2:
        raise ValueError(f"features file must have at least one feature column. Found columns: {list(df.columns)}")
    sample_ids = df.iloc[:, 0].astype(str)
    features = df.iloc[:, 1:].to_numpy(dtype=float)
    return sample_ids, features


def main():
    features_path = "$features"
    algorithm = "$algorithm"
    n_clusters = int("$n_clusters")
    dbscan_eps = float("$dbscan_eps")
    dbscan_min_samples = int("$dbscan_min_samples")
    prefix = "${task.ext.prefix ?: meta.id}"

    sample_ids, x = load_features(features_path)

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
