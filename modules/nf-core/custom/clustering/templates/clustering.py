#!/usr/bin/env python3

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
import yaml
from sklearn.cluster import DBSCAN, KMeans


def parse_eigenvec(path):
    """Parse a PLINK2 .eigenvec file into (sample_ids: pd.Series, pcs: np.ndarray).

    Handles both FID/IID and IID-only header layouts and tolerates a leading
    '#' on the header line that PLINK2 writes by default.
    """
    rows = []
    n_pcs = 0
    mode = None
    with Path(path).open("r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if parts[0].startswith("#"):
                header = [p.lstrip("#") for p in parts]
                if len(header) >= 2 and header[0].upper() == "FID" and header[1].upper() == "IID":
                    mode = "fid_iid"
                elif header[0].upper() == "IID":
                    mode = "iid_only"
                continue
            if mode is None:
                # No header line; infer from the first data row.
                try:
                    float(parts[1])
                    mode = "iid_only"
                except (ValueError, IndexError):
                    mode = "fid_iid"
            if mode == "fid_iid":
                if len(parts) < 3:
                    continue
                rows.append((parts[1], parts[2:]))
            else:
                if len(parts) < 2:
                    continue
                rows.append((parts[0], parts[1:]))
            if n_pcs == 0:
                n_pcs = len(rows[-1][1])
    if not rows:
        raise ValueError(f"No data rows found in eigenvec file {path}")
    sample_ids = pd.Series([r[0] for r in rows], dtype=str)
    pcs = np.array([[float(v) for v in r[1]] for r in rows])
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
            "n_clusters_found": int(len(set(labels)) - (1 if -1 in labels else 0)),
            "n_noise": int(np.sum(labels == -1)),
        }
    else:
        raise ValueError(f"Unknown algorithm '{algorithm}' (expected 'kmeans' or 'dbscan')")

    info["n_samples"] = int(x.shape[0])
    info["n_features"] = int(x.shape[1])

    pd.DataFrame({"sample_id": sample_ids, "cluster": labels}).to_csv(f"{prefix}_clusters.csv", index=False)
    Path(f"{prefix}_clustering_info.json").write_text(json.dumps(info, indent=2))

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
