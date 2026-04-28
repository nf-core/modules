#!/usr/bin/env python3
"""
Unified PLINK2 PCA + Clustering pipeline.

This script:
  1. Converts PLINK2 .eigenvec file to a clean TSV format
  2. Performs clustering (KMeans or DBSCAN) on the principal components
  3. Generates cluster assignments and metadata

Author: Donald Baku (athor)
Date: April 2026
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN

PC_COL_RE = re.compile(r"^PC\d+$", re.IGNORECASE)


def convert_eigenvec_to_tsv(eigenvec_path: Path, out_pca: Path, id_mode: str = "iid"):
    """Convert PLINK2 .eigenvec file to clean TSV (sample_id + PC1 + PC2 + ...)."""
    rows = []
    n_pcs = 0
    mode = None  # "fid_iid" or "iid_only"

    with eigenvec_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()

            # Detect and skip header
            if parts[0].startswith("#"):
                header = [p.lstrip("#") for p in parts]
                if len(header) >= 2 and header[0].upper() == "FID" and header[1].upper() == "IID":
                    mode = "fid_iid"
                elif header[0].upper() == "IID":
                    mode = "iid_only"
                continue

            # Infer mode if there is no header
            if mode is None:
                # Heuristic: if second field is numeric, assume IID-only format
                try:
                    float(parts[1])
                    mode = "iid_only"
                except (ValueError, IndexError):
                    mode = "fid_iid"

            if mode == "fid_iid":
                if len(parts) < 3:
                    continue
                fid = parts[0]
                iid = parts[1]
                pcs = parts[2:]
                sample_id = iid if id_mode == "iid" else f"{fid}:{iid}"

            elif mode == "iid_only":
                if len(parts) < 2:
                    continue
                iid = parts[0]
                pcs = parts[1:]
                sample_id = iid

            else:
                raise ValueError(f"Unrecognized eigenvec format in {eigenvec_path}")

            if n_pcs == 0:
                n_pcs = len(pcs)

            rows.append((sample_id, pcs))

    if not rows:
        raise ValueError(f"No valid data found in {eigenvec_path}")

    header = ["sample_id"] + [f"PC{i+1}" for i in range(n_pcs)]
    with out_pca.open("w") as f:
        f.write("\t".join(header) + "\n")
        for sample_id, pcs in rows:
            f.write(sample_id + "\t" + "\t".join(pcs) + "\n")

    print(f"[INFO] Converted {len(rows)} samples with {n_pcs} PCs → {out_pca}")
    return n_pcs


def read_table_robust(path: str) -> pd.DataFrame:
    """Read TSV file with robust handling of duplicate headers."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    print(f"[DEBUG] Initial read: {df.shape[0]} rows × {df.shape[1]} cols", flush=True)

    col_names_upper = set(str(c).upper() for c in df.columns)

    def is_header_row(row) -> bool:
        row_values_upper = [str(v).upper() for v in row.values]
        overlap = sum(1 for v in row_values_upper if v in col_names_upper)
        if overlap >= 3:
            return True
        header_keywords = {"FID", "IID", "PC1", "PC2", "PC3"}
        if sum(1 for v in row_values_upper if v in header_keywords) >= 2:
            return True
        return False

    bad_rows = df.apply(is_header_row, axis=1)
    if bad_rows.any():
        n_bad = int(bad_rows.sum())
        print(f"[INFO] Removed {n_bad} duplicate header row(s)", flush=True)
        df = df[~bad_rows].copy().reset_index(drop=True)

    print(f"[INFO] After cleanup: {df.shape[0]} rows × {df.shape[1]} cols", flush=True)
    return df


def build_sample_id(df: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """Extract sample IDs and separate numeric feature columns."""

    cols = list(df.columns)

    # 1. Best case: already normalised
    if "sample_id" in df.columns:
        return df["sample_id"].astype(str), df.drop(columns=["sample_id"])

    # 2. PLINK-style IID column
    iid_candidates = [c for c in cols if str(c).upper() == "IID"]
    if iid_candidates:
        iid = iid_candidates[0]
        return df[iid].astype(str), df.drop(columns=[iid])

    # 3. FID + IID case
    fid_candidates = [c for c in cols if str(c).upper() == "FID"]
    if fid_candidates and iid_candidates:
        fid = fid_candidates[0]
        iid = iid_candidates[0]
        sample_ids = df[iid].astype(str)
        return sample_ids, df.drop(columns=[c for c in [fid, iid] if c in df.columns])

    # 4. Fallback: first non-PC column
    pc_cols = [c for c in cols if PC_COL_RE.match(str(c))]
    non_pc_cols = [c for c in cols if c not in pc_cols]
    if non_pc_cols:
        id_col = non_pc_cols[0]
        return df[id_col].astype(str), df.drop(columns=[id_col])

    # 5. Last resort
    return pd.Series([f"sample_{i}" for i in range(len(df))], index=df.index), df

def _is_nullish(x) -> bool:
    if x is None:
        return True
    s = str(x).strip().lower()
    return s in {"", "null", "none", "nan"}


def main():
    ap = argparse.ArgumentParser(description="PLINK2 PCA + Clustering (KMeans / DBSCAN)")
    ap.add_argument("--eigenvec", required=True, help="Input .eigenvec file from PLINK2 --pca")
    ap.add_argument("--algorithm", choices=["kmeans", "dbscan"], default="kmeans")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--dbscan-eps", type=float, default=0.5)
    ap.add_argument("--dbscan-min-samples", type=int, default=5)
    ap.add_argument("--n_init", type=int, default=100)
    ap.add_argument("--init-method", choices=["k-means++", "random"], default="random")
    ap.add_argument("--out-prefix", default="clustering")
    ap.add_argument("--id-mode", choices=["fid_iid", "iid"], default="iid")

    args = ap.parse_args()

    prefix = args.out_prefix

    # Step 1: Convert PLINK2 .eigenvec → clean TSV
    pca_tsv = Path(f"{prefix}_pca_scores.tsv")
    convert_eigenvec_to_tsv(Path(args.eigenvec), pca_tsv, args.id_mode)

    # Step 2: Clustering
    print(f"\n{'='*60}")
    print("CLUSTERING PIPELINE")
    print(f"{'='*60}\n")

    df = read_table_robust(str(pca_tsv))
    sample_ids, df_feats = build_sample_id(df)

    pc_cols = [c for c in df_feats.columns if PC_COL_RE.match(str(c))]
    if not pc_cols:
        raise ValueError("No PC columns found in input")

    X = df_feats[pc_cols].apply(pd.to_numeric, errors="coerce").values

    if np.isnan(X).any():
        raise ValueError("NaN values detected in PCA data")

    print(f"[INFO] Loaded {X.shape[0]} samples × {X.shape[1]} principal components", flush=True)

    # Run clustering
    if args.algorithm == "kmeans":
        print(f"[INFO] Running K-Means (k={args.k})...")
        model = KMeans(
            n_clusters=args.k,
            init=args.init_method,
            n_init=args.n_init,
            random_state=42
        )
        labels = model.fit_predict(X)
        info = {
            "algorithm": "kmeans",
            "k": int(args.k),
            "inertia": float(model.inertia_)
        }
    else:
        print(f"[INFO] Running DBSCAN (eps={args.dbscan_eps}, min_samples={args.dbscan_min_samples})...")
        model = DBSCAN(eps=args.dbscan_eps, min_samples=args.dbscan_min_samples)
        labels = model.fit_predict(X)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = int(np.sum(labels == -1))
        info = {
            "algorithm": "dbscan",
            "eps": float(args.dbscan_eps),
            "min_samples": int(args.dbscan_min_samples),
            "n_clusters_found": int(n_clusters),
            "n_noise": n_noise
        }

    # Save outputs
    out_clusters = f"{prefix}_clusters.csv"
    out_info = f"{prefix}_clustering_info.json"

    pd.DataFrame({"sample_id": sample_ids.astype(str), "cluster": labels}).to_csv(out_clusters, index=False)

    info.update({
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "feature_names": pc_cols,
        "input_file": Path(args.eigenvec).name
    })

    Path(out_info).write_text(json.dumps(info, indent=2))

    print(f"\n[SUCCESS] Clustering completed successfully:")
    print(f"   → Clusters : {out_clusters}")
    print(f"   → Info     : {out_info}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
