#!/usr/bin/env python3

import json
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)


def load_features(path):
    """Read a TSV of `sample_id` + numeric feature columns, indexed by sample_id."""
    df = pd.read_csv(path, sep="\\t")
    if "sample_id" not in df.columns:
        raise ValueError(f"features file must have a 'sample_id' column. Found: {list(df.columns)}")
    df["sample_id"] = df["sample_id"].astype(str)
    return df.set_index("sample_id").apply(pd.to_numeric, errors="coerce").fillna(0.0)


def load_clusters(path):
    """Read a CSV of `sample_id` + `cluster`, returning a Series of int labels."""
    df = pd.read_csv(path)
    if "sample_id" not in df.columns or "cluster" not in df.columns:
        raise ValueError(f"clusters file must have 'sample_id' and 'cluster' columns. Found: {list(df.columns)}")
    df["sample_id"] = df["sample_id"].astype(str)
    return df.set_index("sample_id")["cluster"].astype(int)


def cluster_quality(x, labels):
    """Silhouette / Calinski-Harabasz / Davies-Bouldin for given (x, labels).

    Treats label -1 as DBSCAN noise and excludes those points. Returns None
    for each score when fewer than 2 clusters of more than one point remain.
    """
    mask = labels != -1
    x, labels = x[mask], labels[mask]
    n = len(set(labels))
    valid = 2 <= n < len(x)
    return {
        "silhouette": float(silhouette_score(x, labels)) if valid else None,
        "calinski_harabasz": float(calinski_harabasz_score(x, labels)) if valid else None,
        "davies_bouldin": float(davies_bouldin_score(x, labels)) if valid else None,
    }


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

    joined = load_features(features).join(load_clusters(clusters_path), how="inner")
    if len(joined) < 2:
        raise ValueError(f"Need at least 2 samples with matching sample_id in both inputs. Got {len(joined)}.")

    labels = joined["cluster"].values
    x = joined.drop(columns=["cluster"]).to_numpy(dtype=float)

    # Quality metrics on the supplied labels.
    selected = {"n_clusters": len(set(labels) - {-1}), **cluster_quality(x, labels)}
    pd.DataFrame([selected]).to_csv(f"{prefix}_metrics.tsv", sep="\\t", index=False)
    with open(f"{prefix}_selected.json", "w") as fh:
        json.dump(selected, fh, indent=2)

    # KMeans k-sweep for downstream comparison.
    rows = []
    for k in range(k_min, min(k_max, len(x)) + 1):
        model = KMeans(n_clusters=k, n_init=10, random_state=42).fit(x)
        rows.append({"k": k, "inertia": float(model.inertia_), **cluster_quality(x, model.labels_)})

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
    with open("versions.yml", "w") as fh:
        yaml.dump(versions, fh, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    main()
