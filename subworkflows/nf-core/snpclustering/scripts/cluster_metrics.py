#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score


def _normalise_id_column(df: pd.DataFrame) -> pd.DataFrame:
    cols_upper = {c.upper(): c for c in df.columns}

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
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = _normalise_id_column(df)
    sample_ids = df["sample_id"].astype(str)
    X = df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce")
    X = X.fillna(X.mean(axis=0))
    return X, sample_ids


def load_clusters(path: str) -> pd.Series:
    df = pd.read_csv(path, sep=",", dtype=str)
    df = _normalise_id_column(df)
    if "cluster" not in df.columns:
        raise ValueError("clusters CSV must have a 'cluster' column")
    return df.set_index(df["sample_id"].astype(str))["cluster"].astype(int)


def safe_cluster_metrics(X: np.ndarray, labels: np.ndarray) -> dict:
    uniq = np.unique(labels)
    n_clusters = len(uniq) - (1 if -1 in uniq else 0)
    if n_clusters < 2:
        return {
            "n_clusters": int(n_clusters),
            "silhouette": None,
            "calinski_harabasz": None,
            "davies_bouldin": None,
        }

    mask = labels != -1
    X_use, y_use = X[mask], labels[mask]
    if len(np.unique(y_use)) < 2:
        return {
            "n_clusters": int(n_clusters),
            "silhouette": None,
            "calinski_harabasz": None,
            "davies_bouldin": None,
        }

    return {
        "n_clusters": int(n_clusters),
        "silhouette": float(silhouette_score(X_use, y_use)),
        "calinski_harabasz": float(calinski_harabasz_score(X_use, y_use)),
        "davies_bouldin": float(davies_bouldin_score(X_use, y_use)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--k-min", type=int, default=2)
    ap.add_argument("--k-max", type=int, default=12)
    ap.add_argument("--out-k-sweep", required=True)
    ap.add_argument("--out-selected", required=True)
    ap.add_argument("--out-prefix", required=True)
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
    labels = clusters.loc[common.values].values

    selected = safe_cluster_metrics(X, labels)
    selected["input_clusters"] = Path(args.clusters).name
    selected["input_features"] = Path(args.features).name

    metrics_tsv = f"{args.out_prefix}_metrics.tsv"
    pd.DataFrame([selected]).to_csv(metrics_tsv, sep="\t", index=False)

    kmin, kmax = int(args.k_min), int(args.k_max)
    rows = []
    for k in range(kmin, kmax + 1):
        model = KMeans(n_clusters=k, n_init="auto", random_state=42)
        y = model.fit_predict(X)
        sil = ch = db = None
        if k >= 2:
            sil = float(silhouette_score(X, y))
            ch = float(calinski_harabasz_score(X, y))
            db = float(davies_bouldin_score(X, y))
        rows.append({
            "k": k,
            "inertia": float(model.inertia_),
            "silhouette": sil,
            "calinski_harabasz": ch,
            "davies_bouldin": db,
        })

    sweep_df = pd.DataFrame(rows)
    sweep_df.to_csv(args.out_k_sweep, sep=",", index=False)
    Path(args.out_selected).write_text(json.dumps(selected, indent=2))

    pfx = args.out_prefix
    try:
        import matplotlib.pyplot as plt

        def plot_curve(metric, title, ylabel, out_png):
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

        plot_curve("inertia", "Elbow method (KMeans inertia)", "inertia", f"{pfx}_elbow.png")
        plot_curve("silhouette", "Silhouette score (higher is better)", "silhouette", f"{pfx}_silhouette.png")
        plot_curve("davies_bouldin", "Davies-Bouldin index (lower is better)", "davies_bouldin", f"{pfx}_davies_bouldin.png")
        plot_curve("calinski_harabasz", "Calinski-Harabasz index (higher is better)", "calinski_harabasz", f"{pfx}_calinski.png")

    except Exception as e:
        Path("plot_warning.txt").write_text(f"Plotting failed: {e}\n")


if __name__ == "__main__":
    main()