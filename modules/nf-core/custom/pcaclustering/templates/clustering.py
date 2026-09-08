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
    """Read a sample-by-feature matrix.

    Accepts a generic TSV/TXT (first column = sample IDs, remaining columns
    numeric features) or a PLINK2 ``.eigenvec`` file (``#FID IID PC1 ...``).
    PLINK family/source ID columns are dropped and ``IID`` is used as
    ``sample_id``, so this module can consume ``plink2/pca`` output directly.

    Returns (sample_ids: pd.Series, features: np.ndarray).
    """
    df = pd.read_csv(path, sep=r"\\s+", engine="python", dtype=str)
    if df.empty or df.shape[1] < 2:
        raise ValueError(
            f"features file must have an ID column and at least one feature column. Found columns: {list(df.columns)}"
        )

    df.columns = [str(col).lstrip("#") for col in df.columns]

ignore_cols = {"FID", "IID", "SID", "sample_id"}
    if "IID" in df.columns:
        sample_ids = df["IID"]
        feature_cols = [col for col in df.columns if col not in plink_id_cols]
    elif "sample_id" in df.columns:
        sample_ids = df["sample_id"]
        feature_cols = [col for col in df.columns if col != "sample_id"]
    else:
        sample_ids = df.iloc[:, 0]
        feature_cols = list(df.columns[1:])

    if not feature_cols:
        raise ValueError(
            f"no numeric feature columns left after dropping ID columns. Found columns: {list(df.columns)}"
        )

    features = df[feature_cols].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    return sample_ids.astype(str), features


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
