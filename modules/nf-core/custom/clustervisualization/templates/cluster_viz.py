#!/usr/bin/env python3

import argparse
import os
import shlex

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
import yaml
from sklearn.manifold import TSNE


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


def embed(x, method, umap_neighbors, tsne_perplexity):
    """Project x to 2D using UMAP or t-SNE.

    n_neighbors / perplexity are clamped against sample count so tiny test
    inputs (where the user-specified defaults exceed n_samples) still run.
    """
    n = len(x)
    if method == "umap":
        reducer = umap.UMAP(n_components=2, n_neighbors=min(umap_neighbors, max(2, n - 1)), random_state=42)
    elif method == "tsne":
        reducer = TSNE(n_components=2, perplexity=min(tsne_perplexity, max(2, n - 1)), random_state=42)
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


def main():
    features = "$features"
    clusters_path = "$clusters"
    prefix = "${task.ext.prefix ?: meta.id}"

    # Optional configuration via task.ext.args (nf-core convention).
    raw_args = "$task.ext.args"
    parser = argparse.ArgumentParser()
    parser.add_argument("--umap-neighbors", type=int, default=15)
    parser.add_argument("--tsne-perplexity", type=int, default=30)
    opts = parser.parse_args(shlex.split(raw_args) if raw_args and raw_args != "null" else [])

    joined = load_features(features).join(load_clusters(clusters_path), how="inner")
    if len(joined) < 2:
        raise ValueError(f"Need at least 2 samples with matching sample_id in both inputs. Got {len(joined)}.")

    labels = joined["cluster"].values
    x = joined.drop(columns=["cluster"]).to_numpy(dtype=float)
    sample_ids = joined.index.to_numpy()

    for method in ("umap", "tsne"):
        emb = embed(x, method, opts.umap_neighbors, opts.tsne_perplexity)
        pd.DataFrame({"sample_id": sample_ids, "Dim1": emb[:, 0], "Dim2": emb[:, 1], "cluster": labels}).to_csv(
            f"{prefix}.{method}.tsv", sep="\\t", index=False
        )
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
    with open("versions.yml", "w") as fh:
        yaml.dump(versions, fh, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    main()
