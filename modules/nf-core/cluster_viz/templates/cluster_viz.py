#!/usr/bin/env python3

"""Cluster visualizations.

Produces three 2D plots, all colored by cluster label:
  - PCA (first two columns from pca_scores)
  - UMAP (computed on the feature matrix used for clustering)
  - t-SNE (computed on the feature matrix used for clustering)

Also writes UMAP and t-SNE coordinates to TSV.
"""

import argparse

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE


def _normalise_id_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Handles the header formats that FlashPCA/PLINK2 produces:
      - '#IID' (PLINK2 eigenvec: leading hash on first column)
      - 'IID'  (FlashPCA / older PLINK)
      - 'FID', 'IID' (two-column prefix)
      - 'sample_id' (already normalised)
    """
    # Strip leading '#' (PLINK2 eigenvec writes '#IID' as the first column)
    df = df.rename(columns=lambda c: c.lstrip("#"))

    cols_upper = {c.upper(): c for c in df.columns}

    # Remove duplicate header row (IID value == "FID" or "IID")
    if "IID" in cols_upper:
        iid_col = cols_upper["IID"]
        dup_mask = df[iid_col].str.upper().isin({"FID", "IID"})
        if dup_mask.any():
            df = df[~dup_mask].copy().reset_index(drop=True)

    cols_upper = {c.upper(): c for c in df.columns}

    if "SAMPLE_ID" in cols_upper:
        return df

    if "IID" in cols_upper:
        iid_col = cols_upper["IID"]
        iid_numeric = pd.to_numeric(df[iid_col], errors="coerce").notna().all()
        if iid_numeric:
            df = df.drop(columns=[iid_col])
            df = df.rename(columns={df.columns[0]: "sample_id"})
        else:
            df = df.rename(columns={iid_col: "sample_id"})

        fid_cols = [c for c in df.columns if c.upper() == "FID"]
        if fid_cols:
            df = df.drop(columns=fid_cols)
        return df

    raise ValueError(
        f"Cannot find sample ID column (expected 'sample_id' or 'IID'). "
        f"Found: {list(df.columns)}"
    )


def load_features(path: str) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
    df = _normalise_id_column(df)
    sample_ids = df["sample_id"].astype(str)
    X = (
        df.drop(columns=["sample_id"])
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
    )
    return X, sample_ids


def load_clusters(path: str) -> pd.Series:
    df = pd.read_csv(path, sep=",", dtype=str)
    df = _normalise_id_column(df)
    if "cluster" not in df.columns:
        raise ValueError("clusters CSV must have a 'cluster' column")
    return df.set_index("sample_id")["cluster"].astype(int)


def safe_perplexity(n_samples: int, requested: float) -> float:
    if n_samples <= 3:
        return 1.0
    upper = (n_samples - 1) / 3.0
    return float(max(2.0, min(requested, upper)))


def compute_umap(X: np.ndarray, n_neighbors: int, min_dist: float) -> np.ndarray:
    try:
        import umap
        return umap.UMAP(
            n_components=2,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=42,
        ).fit_transform(X)
    except Exception as e:
        print(f"[WARN] UMAP failed, fallback to first 2 feature columns: {e}")
        if X.shape[1] >= 2:
            return X[:, :2]
        elif X.shape[1] == 1:
            return np.column_stack([X[:, 0], np.zeros(X.shape[0])])
        else:
            return np.zeros((X.shape[0], 2))


def compute_tsne(X: np.ndarray, perplexity: float, max_iter: int) -> np.ndarray:
    return TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        random_state=42,
        max_iter=max_iter,
        learning_rate="auto",
    ).fit_transform(X)


def plot_scatter(df, x, y, out_png, title, xlabel=None, ylabel=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.figure(figsize=(7, 5))
    labels = df["cluster"].astype(int).values
    uniq = np.unique(labels)
    sc = plt.scatter(df[x], df[y], c=labels, cmap="Paired", s=24,linewidths=0.4, alpha=0.85)
    plt.title(title)
    plt.xlabel(xlabel or x)
    plt.ylabel(ylabel or y)
    plt.grid(True, alpha=0.5)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=7,
               markerfacecolor=sc.cmap(sc.norm(k)), markeredgecolor="none",
               label=f"Cluster {k}")
        for k in uniq
    ]
    plt.legend(handles=handles, title="Clusters", loc="center left",
               bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0, frameon=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser(description="PCA + UMAP + t-SNE plots colored by cluster")
    ap.add_argument("--features", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--pca-scores", required=True)
    ap.add_argument("--tsne-perplexity", type=float, default=30.0)
    ap.add_argument("--tsne-iter", type=int, default=1000)
    ap.add_argument("--umap-neighbors", type=int, default=15)
    ap.add_argument("--umap-min-dist", type=float, default=0.1)
    ap.add_argument("--out-umap-tsv", required=True)
    ap.add_argument("--out-tsne-tsv", required=True)
    ap.add_argument("--out-umap-png", required=True)
    ap.add_argument("--out-tsne-png", required=True)
    ap.add_argument("--out-pca-png", required=True)
    args = ap.parse_args()

    X_df, sample_ids = load_features(args.features)
    clusters = load_clusters(args.clusters)

    common = sample_ids[sample_ids.isin(clusters.index)]
    if len(common) == 0:
        raise ValueError(
            f"No overlapping sample_id between features and clusters.\n"
            f"  features IDs (first 5): {sample_ids.head().tolist()}\n"
            f"  clusters IDs (first 5): {list(clusters.index[:5])}"
        )

    X = X_df.loc[common.index].values
    y = clusters.loc[common.values].values

    umap_coords = compute_umap(X, args.umap_neighbors, args.umap_min_dist)
    umap_df = pd.DataFrame({
        "sample_id": common.values,
        "x": umap_coords[:, 0],
        "y": umap_coords[:, 1],
        "cluster": y,
    })
    umap_df.to_csv(args.out_umap_tsv, sep="\t", index=False)
    plot_scatter(umap_df, "x", "y", args.out_umap_png, "UMAP embedding")

    perp = safe_perplexity(len(common), args.tsne_perplexity)
    tsne_coords = compute_tsne(X, perp, args.tsne_iter)
    tsne_df = pd.DataFrame({
        "sample_id": common.values,
        "x": tsne_coords[:, 0],
        "y": tsne_coords[:, 1],
        "cluster": y,
    })
    tsne_df.to_csv(args.out_tsne_tsv, sep="\t", index=False)
    plot_scatter(tsne_df, "x", "y", args.out_tsne_png,
                 f"t-SNE (perplexity={perp:.1f})")

    pca_df = pd.read_csv(args.pca_scores, sep=r"\s+", engine="python", dtype=str)
    pca_df = _normalise_id_column(pca_df)
    comp_cols = [c for c in pca_df.columns if c != "sample_id"]
    if len(comp_cols) < 2:
        raise ValueError("pca_scores must have at least 2 PC columns")
    c1, c2 = comp_cols[0], comp_cols[1]
    for col in [c1, c2]:
        pca_df[col] = pd.to_numeric(pca_df[col], errors="coerce")
    merged = pca_df.merge(umap_df[["sample_id", "cluster"]], on="sample_id", how="inner")
    plot_scatter(merged, c1, c2, args.out_pca_png,
                 "PCA", c1, c2)


if __name__ == "__main__":
    main()
