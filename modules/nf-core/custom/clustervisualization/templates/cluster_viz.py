#!/usr/bin/env python3

import os

# numba (UMAP) and matplotlib write caches; redirect to /tmp so the script works
# inside read-only container filesystems.
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp")
os.environ.setdefault("MPLCONFIGDIR", "/tmp")

import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn
import umap
from sklearn.manifold import TSNE


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


def embed(x, method):
    if method == "umap":
        # min(15, n-1) keeps UMAP working on tiny test inputs where the default
        # n_neighbors (15) exceeds sample count.
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=min(15, max(2, len(x) - 1)),
            random_state=42,
        )
    elif method == "tsne":
        reducer = TSNE(
            n_components=2,
            perplexity=min(30, max(2, len(x) - 1)),
            random_state=42,
        )
    else:
        raise ValueError(f"Unknown method '{method}' (expected 'umap' or 'tsne')")
    return reducer.fit_transform(x)


def plot_embedding(emb, labels, method, out_png):
    plt.figure(figsize=(8, 6))
    palette = sns.color_palette("tab10", n_colors=max(1, len(np.unique(labels))))
    sns.scatterplot(
        x=emb[:, 0],
        y=emb[:, 1],
        hue=labels.astype(str),
        palette=palette,
        alpha=0.8,
        s=60,
        edgecolor="k",
        linewidth=0.3,
    )
    plt.title(f"{method.upper()} projection colored by cluster")
    plt.xlabel(f"{method.upper()} 1")
    plt.ylabel(f"{method.upper()} 2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def write_embedding(emb, sample_ids, labels, out_tsv):
    df = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "Dim1": emb[:, 0],
            "Dim2": emb[:, 1],
            "cluster": labels,
        }
    )
    df.to_csv(out_tsv, sep="\\t", index=False)


def main():
    features = "$features"
    clusters_path = "$clusters"
    prefix = "${task.ext.prefix ?: meta.id}"

    x, sample_ids = load_features(features)
    clusters = load_clusters(clusters_path)

    common = sample_ids[sample_ids.isin(clusters.index)]
    if len(common) < 2:
        raise ValueError(f"Need at least 2 samples with matching sample_id in both inputs. Got {len(common)}.")

    x_aligned = x[common.index]
    aligned_ids = common.values
    labels = clusters.loc[aligned_ids].values

    for method in ("umap", "tsne"):
        emb = embed(x_aligned, method)
        write_embedding(emb, aligned_ids, labels, f"{prefix}.{method}.tsv")
        plot_embedding(emb, labels, method, f"{prefix}.{method}.png")

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
            "umap-learn": umap.__version__,
            "scikit-learn": sklearn.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


if __name__ == "__main__":
    main()
