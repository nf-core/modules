#!/usr/bin/env python3

import os

# Fix numba + matplotlib in read-only Singularity container
os.environ["NUMBA_CACHE_DIR"] = "/tmp"
os.environ["MPLCONFIGDIR"] = "/tmp"

import platform

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sklearn
import umap as umap_module
from sklearn.manifold import TSNE
from umap import UMAP

matplotlib.use("Agg")


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


def _normalise_id_column(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).lstrip("#") for c in df.columns]

    cols_upper = {str(c).upper(): c for c in df.columns}

    if "IID" in cols_upper:
        iid_col = cols_upper["IID"]
        dup_mask = df[iid_col].astype(str).str.upper().isin({"FID", "IID"})
        if dup_mask.any():
            df = df.loc[~dup_mask].copy().reset_index(drop=True)

    cols_upper = {str(c).upper(): c for c in df.columns}

    if "SAMPLE_ID" in cols_upper:
        sample_col = cols_upper["SAMPLE_ID"]
        if sample_col != "sample_id":
            df = df.rename(columns={sample_col: "sample_id"})
        return df

    if "IID" in cols_upper:
        iid_col = cols_upper["IID"]
        iid_numeric = pd.to_numeric(df[iid_col], errors="coerce").notna().all()

        if iid_numeric:
            df = df.drop(columns=[iid_col])
            if len(df.columns) == 0:
                raise ValueError("Cannot infer sample_id after dropping numeric IID column")
            df = df.rename(columns={df.columns[0]: "sample_id"})
        else:
            df = df.rename(columns={iid_col: "sample_id"})

        fid_cols = [c for c in df.columns if str(c).upper() == "FID"]
        if fid_cols:
            df = df.drop(columns=fid_cols)

        return df

    raise ValueError(f"Cannot find sample ID column (expected 'sample_id' or 'IID'). Found: {list(df.columns)}")


def load_features(path: str) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path, sep="\\t", dtype=str)
    df = _normalise_id_column(df)

    if "sample_id" not in df.columns:
        raise ValueError("features file must contain a sample_id column after normalization")

    sample_ids = df["sample_id"].astype(str)
    x = df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.mean(numeric_only=True))
    x = x.fillna(0.0)

    return x, sample_ids


def load_clusters(path: str) -> tuple[pd.DataFrame, str]:
    """Load clusters and return (df, mode). Same logic as cluster_metrics."""
    df = pd.read_csv(path, sep=",", dtype=str)
    df = df.copy()
    df.columns = [str(c).lstrip("#") for c in df.columns]

    cols_upper = {str(c).upper(): c for c in df.columns}

    if "CLUSTER" not in cols_upper:
        raise ValueError("clusters CSV must have a 'cluster' column")

    cluster_col = cols_upper["CLUSTER"]

    if "SAMPLE_ID" in cols_upper:
        sample_col = cols_upper["SAMPLE_ID"]
        out = df[[sample_col, cluster_col]].copy()
        out.columns = ["sample_id", "cluster"]
        out["sample_id"] = out["sample_id"].astype(str)
        out["cluster"] = pd.to_numeric(out["cluster"], errors="raise").astype(int)
        return out, "sample_id"

    try:
        norm = _normalise_id_column(df.copy())
        if "sample_id" in norm.columns and "cluster" in norm.columns:
            out = norm[["sample_id", "cluster"]].copy()
            out["sample_id"] = out["sample_id"].astype(str)
            out["cluster"] = pd.to_numeric(out["cluster"], errors="raise").astype(int)
            return out, "sample_id"
    except Exception:
        pass

    other_cols = [c for c in df.columns if c != cluster_col]

    if len(other_cols) == 1:
        candidate = other_cols[0]
        candidate_vals = df[candidate].astype(str)

        if not (
            len(candidate_vals) > 0 and float(pd.to_numeric(candidate_vals, errors="coerce").notna().mean()) >= 0.8
        ):
            out = pd.DataFrame(
                {
                    "sample_id": candidate_vals,
                    "cluster": pd.to_numeric(df[cluster_col], errors="raise").astype(int),
                }
            )
            return out, "sample_id"

    out = pd.DataFrame({"cluster": pd.to_numeric(df[cluster_col], errors="raise").astype(int)})
    return out, "row_order"


def plot_embedding(x: np.ndarray, labels: np.ndarray, method: str, prefix: str) -> None:
    """Plot UMAP or t-SNE with cluster coloring."""
    if method == "umap":
        reducer = UMAP(random_state=42)
        embedding = reducer.fit_transform(x)
        title = "UMAP"
        out_tsv = f"{prefix}.umap.tsv"
        out_png = f"{prefix}.umap.png"
    else:  # tsne
        reducer = TSNE(n_components=2, random_state=42, perplexity=min(30, len(x) - 1))
        embedding = reducer.fit_transform(x)
        title = "t-SNE"
        out_tsv = f"{prefix}.tsne.tsv"
        out_png = f"{prefix}.tsne.png"

    # Save embedding
    emb_df = pd.DataFrame(embedding, columns=["Dim1", "Dim2"])
    emb_df["cluster"] = labels
    emb_df.to_csv(out_tsv, sep="\\t", index=False)

    # Plot
    plt.figure(figsize=(8, 6))
    palette = sns.color_palette("tab10", n_colors=len(np.unique(labels)))
    sns.scatterplot(
        x=embedding[:, 0],
        y=embedding[:, 1],
        hue=labels.astype(str),
        palette=palette,
        alpha=0.8,
        s=60,
        edgecolor="k",
        linewidth=0.3,
    )
    plt.title(f"{title} projection of features colored by cluster")
    plt.xlabel(f"{title} 1")
    plt.ylabel(f"{title} 2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def main() -> None:
    features = "$features"
    clusters_path = "$clusters"
    prefix = "${task.ext.prefix ?: meta.id}"

    x_df, sample_ids = load_features(features)
    clusters_df, cluster_mode = load_clusters(clusters_path)

    if cluster_mode == "sample_id":
        clusters = clusters_df.set_index("sample_id")["cluster"]
        common = sample_ids[sample_ids.isin(clusters.index)]
        if len(common) > 0:
            x = x_df.loc[common.index].values
            labels = clusters.loc[common.values].values
        elif len(clusters_df) == len(sample_ids):
            x = x_df.values
            labels = clusters_df["cluster"].values
        else:
            raise ValueError("No overlapping sample_id between features and clusters")
    else:
        if len(clusters_df) != len(sample_ids):
            raise ValueError("Row counts do not match and no sample_id column found")
        x = x_df.values
        labels = clusters_df["cluster"].values

    if len(x) < 2:
        raise ValueError("Need at least 2 samples for embedding")

    # Generate both embeddings
    plot_embedding(x, labels, "umap", prefix)
    plot_embedding(x, labels, "tsne", prefix)

    # versions.yml
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
            "umap-learn": umap_module.__version__,
            "scikit-learn": sklearn.__version__,
        }
    }
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))


if __name__ == "__main__":
    main()
