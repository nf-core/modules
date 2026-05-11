#!/usr/bin/env python3

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import sklearn
from sklearn.cluster import KMeans
from sklearn.metrics import (
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)

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
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = _normalise_id_column(df)

    if "sample_id" not in df.columns:
        raise ValueError("features file must contain a sample_id column after normalization")

    sample_ids = df["sample_id"].astype(str)
    x = df.drop(columns=["sample_id"]).apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.mean(numeric_only=True))
    x = x.fillna(0.0)

    return x, sample_ids


def _looks_mostly_numeric(s: pd.Series) -> bool:
    if len(s) == 0:
        return False
    parsed = pd.to_numeric(s.astype(str), errors="coerce")
    return float(parsed.notna().mean()) >= 0.8


def load_clusters(path: str) -> tuple[pd.DataFrame, str]:
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

        if not _looks_mostly_numeric(candidate_vals):
            out = pd.DataFrame(
                {
                    "sample_id": candidate_vals,
                    "cluster": pd.to_numeric(df[cluster_col], errors="raise").astype(int),
                }
            )
            return out, "sample_id"

    out = pd.DataFrame({"cluster": pd.to_numeric(df[cluster_col], errors="raise").astype(int)})
    return out, "row_order"


def safe_cluster_metrics(x: np.ndarray, labels: np.ndarray) -> dict:
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
    x_use, y_use = x[mask], labels[mask]

    if len(x_use) < 2 or len(np.unique(y_use)) < 2:
        return {
            "n_clusters": int(n_clusters),
            "silhouette": None,
            "calinski_harabasz": None,
            "davies_bouldin": None,
        }

    return {
        "n_clusters": int(n_clusters),
        "silhouette": float(silhouette_score(x_use, y_use)),
        "calinski_harabasz": float(calinski_harabasz_score(x_use, y_use)),
        "davies_bouldin": float(davies_bouldin_score(x_use, y_use)),
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

    x_df, sample_ids = load_features(args.features)
    clusters_df, cluster_mode = load_clusters(args.clusters)

    if cluster_mode == "sample_id":
        clusters = clusters_df.set_index("sample_id")["cluster"]
        common = sample_ids[sample_ids.isin(clusters.index)]

        if len(common) > 0:
            x = x_df.loc[common.index].values
            labels = clusters.loc[common.values].values
            aligned_ids = common.astype(str).tolist()
            alignment_mode = "sample_id"
        elif len(clusters_df) == len(sample_ids):
            x = x_df.values
            labels = clusters_df["cluster"].values
            aligned_ids = sample_ids.astype(str).tolist()
            alignment_mode = "row_order_fallback"
        else:
            raise ValueError(
                f"No overlapping sample_id between features and clusters.\\n"
                f"  features IDs (first 5): {sample_ids.head().tolist()}\\n"
                f"  clusters IDs (first 5): {list(clusters.index[:5])}"
            )
    else:
        if len(clusters_df) != len(sample_ids):
            raise ValueError(
                "clusters CSV has no usable sample_id column and row counts do not match.\\n"
                f"  n_features={len(sample_ids)}\\n"
                f"  n_clusters={len(clusters_df)}"
            )
        x = x_df.values
        labels = clusters_df["cluster"].values
        aligned_ids = sample_ids.astype(str).tolist()
        alignment_mode = "row_order"

    if len(x) < 2:
        raise ValueError("Need at least 2 samples to compute cluster metrics")

    selected = safe_cluster_metrics(x, labels)
    selected["input_clusters"] = Path(args.clusters).name
    selected["input_features"] = Path(args.features).name
    selected["n_samples_used"] = int(len(aligned_ids))
    selected["alignment_mode"] = alignment_mode

    metrics_tsv = f"{args.out_prefix}_metrics.tsv"
    pd.DataFrame([selected]).to_csv(metrics_tsv, sep="\t", index=False)

    rows = []
    max_k = min(int(args.k_max), len(x))
    for k in range(int(args.k_min), max_k + 1):
        model = KMeans(n_clusters=k, n_init="auto", random_state=42)
        y = model.fit_predict(x)

        sil = ch = db = None
        if 1 < len(np.unique(y)) < len(x):
            sil = float(silhouette_score(x, y))
            ch = float(calinski_harabasz_score(x, y))
            db = float(davies_bouldin_score(x, y))

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
    sweep_df.to_csv(args.out_k_sweep, sep=",", index=False, float_format="%.10g")
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

        if not sweep_df.empty:
            plot_curve("inertia", "Elbow method (KMeans inertia)", "inertia", f"{pfx}_elbow.png")
            plot_curve("silhouette", "Silhouette score (higher is better)", "silhouette", f"{pfx}_silhouette.png")
            plot_curve(
                "davies_bouldin",
                "Davies-Bouldin index (lower is better)",
                "davies_bouldin",
                f"{pfx}_davies_bouldin.png",
            )
            plot_curve(
                "calinski_harabasz",
                "Calinski-Harabasz index (higher is better)",
                "calinski_harabasz",
                f"{pfx}_calinski.png",
            )

    except Exception as e:
        Path("plot_warning.txt").write_text("Plotting failed: " + str(e) + "\\n")

    # === VERSIONS.YML (fix review) ===
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
    prefix = "${task.ext.prefix ?: meta.id}"

    sys.argv = [
        "cluster_metrics.py",
        "--features",
        "$features",
        "--clusters",
        "$clusters",
        "--k-min",
        "2",
        "--k-max",
        "12",
        "--out-k-sweep",
        f"{prefix}_k_sweep.csv",
        "--out-selected",
        f"{prefix}_selected.json",
        "--out-prefix",
        prefix,
    ]

    main()
