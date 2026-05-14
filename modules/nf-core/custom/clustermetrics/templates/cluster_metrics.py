#!/usr/bin/env python3

import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string (nf-core standard)."""
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


def load_features(path):
    df = pd.read_csv(path, sep="\\t")
    if "sample_id" not in df.columns:
        raise ValueError(f"features file must have a 'sample_id' column. Found: {list(df.columns)}")
    sample_ids = df["sample_id"].astype(str)
    x = df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce").fillna(0.0).values
    return x, sample_ids


def load_clusters(path):
    df = pd.read_csv(path)
    if "sample_id" not in df.columns or "cluster" not in df.columns:
        raise ValueError(f"clusters file must have 'sample_id' and 'cluster' columns. Found: {list(df.columns)}")
    return df.set_index(df["sample_id"].astype(str))["cluster"].astype(int)


def cluster_quality(x, labels):
    """Silhouette / Calinski-Harabasz / Davies-Bouldin for given (x, labels).

    Treats label -1 as DBSCAN noise and excludes those points.
    Returns None for each score when fewer than 2 clusters remain.
    """
    uniq = np.unique(labels)
    n_clusters = len(uniq) - (1 if -1 in uniq else 0)
    out = {
        "n_clusters": int(n_clusters),
        "silhouette": None,
        "calinski_harabasz": None,
        "davies_bouldin": None,
    }
    if n_clusters < 2:
        return out
    mask = labels != -1
    x_use, y_use = x[mask], labels[mask]
    if len(x_use) < 2 or len(np.unique(y_use)) < 2:
        return out
    out["silhouette"] = float(silhouette_score(x_use, y_use))
    out["calinski_harabasz"] = float(calinski_harabasz_score(x_use, y_use))
    out["davies_bouldin"] = float(davies_bouldin_score(x_use, y_use))
    return out


def plot_curve(sweep_df, metric, title, ylabel, out_png):
    plt.figure(figsize=(7, 4.5))
    vals = sweep_df[metric].dropna()
    ks = sweep_df.loc[vals.index, "k"]
    plt.plot(ks, vals, marker="o")
    plt.xticks(sweep_df["k"].tolist())
    plt.title(title)
    plt.xlabel("k")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main():
    features = "$features"
    clusters_path = "$clusters"
    prefix = "${task.ext.prefix ?: meta.id}"
    k_min, k_max = 2, 12

    x, sample_ids = load_features(features)
    clusters = load_clusters(clusters_path)

    common = sample_ids[sample_ids.isin(clusters.index)]
    if len(common) < 2:
        raise ValueError(f"Need at least 2 samples with matching sample_id in both inputs. Got {len(common)}.")

    x_aligned = x[common.index]
    labels = clusters.loc[common.values].values

    # Per-cluster quality metrics on the supplied labels.
    selected = cluster_quality(x_aligned, labels)
    pd.DataFrame([selected]).to_csv(f"{prefix}_metrics.tsv", sep="\\t", index=False)
    Path(f"{prefix}_selected.json").write_text(json.dumps(selected, indent=2))

    # KMeans k-sweep for downstream comparison.
    rows = []
    max_k = min(k_max, len(x_aligned))
    for k in range(k_min, max_k + 1):
        model = KMeans(n_clusters=k, n_init=10, random_state=42)
        y = model.fit_predict(x_aligned)
        sil = ch = db = None
        if 1 < len(np.unique(y)) < len(x_aligned):
            sil = float(silhouette_score(x_aligned, y))
            ch = float(calinski_harabasz_score(x_aligned, y))
            db = float(davies_bouldin_score(x_aligned, y))
        rows.append(
            {
                "k": k,
                "inertia": float(model.inertia_),
                "silhouette": sil,
                "calinski_harabasz": ch,
                "davies_bouldin": db,
            }
        )

    sweep_df = pd.DataFrame(rows)
    sweep_df.to_csv(f"{prefix}_k_sweep.csv", index=False, float_format="%.10g")

    if not sweep_df.empty:
        plot_curve(sweep_df, "inertia", "Elbow method (KMeans inertia)", "inertia", f"{prefix}_elbow.png")
        plot_curve(
            sweep_df, "silhouette", "Silhouette score (higher is better)", "silhouette", f"{prefix}_silhouette.png"
        )
        plot_curve(
            sweep_df,
            "davies_bouldin",
            "Davies-Bouldin index (lower is better)",
            "davies_bouldin",
            f"{prefix}_davies_bouldin.png",
        )
        plot_curve(
            sweep_df,
            "calinski_harabasz",
            "Calinski-Harabasz index (higher is better)",
            "calinski_harabasz",
            f"{prefix}_calinski.png",
        )

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "scikit-learn": sklearn.__version__,
            "matplotlib": matplotlib.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


if __name__ == "__main__":
    main()
