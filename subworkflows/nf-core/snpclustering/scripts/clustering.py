#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN

PC_COL_RE = re.compile(r"^PC\d+$", re.IGNORECASE)


def read_table_robust(path: str) -> pd.DataFrame:
    """Leggi TSV file gestendo header multipli."""
    df = pd.read_csv(path, sep="\t", dtype=str)

    print(f"[DEBUG] Initial read: {df.shape[0]} rows × {df.shape[1]} cols", flush=True)

    # Rimuovi righe header duplicate
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
        df = df[~bad_rows].copy()
        df.reset_index(drop=True, inplace=True)

    print(f"[INFO] After cleanup: {df.shape[0]} rows × {df.shape[1]} cols", flush=True)
    return df


def build_sample_id(df: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    """Estrai sample IDs e features."""
    cols = list(df.columns)

    # BUG FIX: Only treat a PC column as sample IDs if it is ALSO the very first
    # column in the file AND all other PC columns parse as numeric. The original
    # check fired whenever PC1 happened to contain strings, which caused it to
    # misidentify legitimate PC1 floating-point values as IDs when the column
    # appeared first (e.g. after a malformed merge dropped the IID column).
    # Now we require that a "PC-as-ID" column be non-numeric AND that at least
    # one other numeric PC column exists, making the detection much more robust.
    pc_cols = [c for c in cols if re.match(r"^PC\d+$", str(c), re.IGNORECASE)]
    if pc_cols:
        first_pc = pc_cols[0]
        first_pc_numeric = pd.to_numeric(df[first_pc], errors="coerce")
        remaining_pc_cols = pc_cols[1:]
        has_numeric_pcs = any(
            pd.to_numeric(df[c], errors="coerce").notna().any()
            for c in remaining_pc_cols
        )
        # Only treat as IDs if: all values non-numeric AND other numeric PCs exist
        if first_pc_numeric.isna().all() and has_numeric_pcs:
            print(
                f"[INFO] Column '{first_pc}' contains sample IDs, using as sample_id",
                flush=True,
            )
            return df[first_pc].astype(str), df.drop(columns=[first_pc])

    # Caso standard: cerca IID esplicito
    iid_candidates = [c for c in cols if str(c).upper() == "IID"]
    if iid_candidates:
        iid = iid_candidates[0]
        return df[iid].astype(str), df.drop(columns=[iid])

    # Prima colonna come ID
    if len(cols) >= 2:
        return df[cols[0]].astype(str), df.drop(columns=[cols[0]])

    # Fallback
    return pd.Series([f"sample_{i}" for i in range(len(df))], index=df.index), df


def _is_nullish(x) -> bool:
    if x is None:
        return True
    s = str(x).strip().lower()
    return s in {"", "null", "none", "nan"}


def _safe_write(path_str: str, content: str) -> bool:
    """
    Write content to a file whose name may contain characters that are valid on
    Linux but that Python's pathlib / open() handle inconsistently when they look
    like glob patterns (square brackets, spaces, commas).

    We use the low-level open() with the exact string so the OS creates the file
    with that literal name, bypassing any glob expansion.

    Returns True on success, False on failure.
    """
    try:
        with open(path_str, "w", encoding="utf-8") as fh:
            fh.write(content)
        return True
    except Exception as e:
        print(f"[WARN] Could not write '{path_str}': {e}", flush=True)
        return False


def main() -> None:
    ap = argparse.ArgumentParser(description="Cluster PCA features")

    ap.add_argument("--features", required=False)
    ap.add_argument("--algorithm", choices=["kmeans", "dbscan"], default=None)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--dbscan-eps", type=float, default=None)
    ap.add_argument("--dbscan-min-samples", type=int, default=None)
    ap.add_argument("--out-clusters", required=False)
    ap.add_argument("--out-info", required=False)
    ap.add_argument("--n_init", type=int, default=100,
                    help="Number of KMeans initializations (default: 100)")
    ap.add_argument("--init-method", choices=["k-means++", "random"], default="random",
                    help="KMeans initialization method (default: random)")

    # Alias Nextflow
    ap.add_argument("--input", nargs="+", default=None)
    ap.add_argument("--method", dest="algorithm_compat", default=None)
    ap.add_argument("--n_clusters", dest="k_compat", type=int, default=None)
    ap.add_argument("--db_eps", dest="db_eps_compat", type=float, default=None)

    # type=str per accettare anche "null"
    ap.add_argument("--db_min_samples", dest="db_min_samples_compat", type=str, default=None)

    ap.add_argument("--out_prefix", nargs="+", default=None)

    args = ap.parse_args()

    # ---------- features (3° file della tupla input) ----------
    if not args.features and args.input:
        if len(args.input) >= 3:
            args.features = args.input[2]
        else:
            ap.error("--input must include at least 3 files")

    # ---------- algorithm: gestisci "null" ----------
    alg = args.algorithm if not _is_nullish(args.algorithm) else None
    alg_compat = args.algorithm_compat if not _is_nullish(args.algorithm_compat) else None
    args.algorithm = alg or alg_compat or "kmeans"

    # ---------- k ----------
    args.k = args.k if args.k is not None else (args.k_compat if args.k_compat is not None else 3)

    # ---------- dbscan params ----------
    args.dbscan_eps = args.dbscan_eps if args.dbscan_eps is not None else (
        args.db_eps_compat if args.db_eps_compat is not None else 0.5
    )

    if args.dbscan_min_samples is None:
        v = args.db_min_samples_compat
        if _is_nullish(v):
            args.dbscan_min_samples = 5
        else:
            args.dbscan_min_samples = int(v)

    # ---------- out_prefix: conserva RAW e CLEAN ----------
    raw_prefix = args.out_prefix
    if isinstance(raw_prefix, list):
        raw_prefix = " ".join(raw_prefix)
    raw_prefix = str(raw_prefix) if raw_prefix is not None else None

    clean_prefix = re.sub(r"[\[\]\s,]+", "_", str(raw_prefix or "clustering")).strip("_") or "clustering"

    # se non hai passato out-clusters/out-info, usa CLEAN di default
    args.out_clusters = args.out_clusters or f"{clean_prefix}_clusters.csv"
    args.out_info = args.out_info or f"{clean_prefix}_clustering_info.json"

    if not args.features:
        ap.error("Missing --features")

    print(f"\n{'='*60}")
    print("CLUSTERING PIPELINE")
    print(f"{'='*60}\n")

    # Leggi dati
    df = read_table_robust(args.features)
    sample_ids, df_feats = build_sample_id(df)

    # Seleziona colonne PC
    pc_cols = [c for c in df_feats.columns if PC_COL_RE.match(str(c))]
    if not pc_cols:
        raise ValueError("No PC columns found in input")

    print(f"[INFO] Using {len(pc_cols)} PC columns: {pc_cols}", flush=True)

    # Converti a numerico
    X_df = df_feats[pc_cols].apply(pd.to_numeric, errors="coerce")

    # Check NaN
    nan_count = int(X_df.isna().sum().sum())
    if nan_count > 0:
        print(f"[ERROR] Found {nan_count} NaN values after conversion!", flush=True)
        print("[DEBUG] NaN per column:", flush=True)
        for col in pc_cols:
            n = int(X_df[col].isna().sum())
            if n > 0:
                print(f"  {col}: {n} NaN", flush=True)

        rows_with_nan = X_df[X_df.isna().any(axis=1)].head(3)
        print("\n[DEBUG] First 3 rows with NaN:", flush=True)
        print(rows_with_nan)
        print("\n[DEBUG] Corresponding raw values:", flush=True)
        print(df_feats.loc[rows_with_nan.index, pc_cols].head(3))

        raise ValueError(f"Cannot proceed: {nan_count} NaN values in data")

    X = X_df.values
    print(f"[INFO] Loaded {X.shape[0]} samples × {X.shape[1]} features (no NaN)\n", flush=True)

    # Clustering
    info = {
        "algorithm": args.algorithm,
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "feature_names": pc_cols,
    }

    if args.algorithm == "kmeans":
        print(f"[INFO] Running K-Means (k={args.k})...", flush=True)
        model = KMeans(n_clusters=args.k, init=args.init_method, n_init=args.n_init, random_state=42)
        labels = model.fit_predict(X)
        info.update({"k": int(args.k), "inertia": float(model.inertia_)})
    else:
        print(f"[INFO] Running DBSCAN (eps={args.dbscan_eps}, min_samples={args.dbscan_min_samples})...", flush=True)
        model = DBSCAN(eps=args.dbscan_eps, min_samples=args.dbscan_min_samples)
        labels = model.fit_predict(X)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = int(np.sum(labels == -1))
        info.update({
            "eps": float(args.dbscan_eps),
            "min_samples": int(args.dbscan_min_samples),
            "n_clusters_found": int(n_clusters),
            "n_noise": n_noise,
        })
        print(f"[INFO] Found {n_clusters} clusters, {n_noise} noise points", flush=True)

    # Salva risultati (CLEAN)
    out_df = pd.DataFrame({"sample_id": sample_ids.astype(str), "cluster": labels})
    out_df.to_csv(args.out_clusters, index=False)
    info_text = json.dumps(info, indent=2) + "\n"
    Path(args.out_info).write_text(info_text)

    # ---- compat Nextflow: se Nextflow si aspetta il RAW filename con parentesi/spazi, scrivilo anche ----
    # BUG FIX: Path() silently mangled bracket filenames due to glob interpretation.
    # We now use _safe_write() which calls open() with the raw string directly,
    # guaranteeing the OS creates the exact literal filename Nextflow expects.
    if raw_prefix and ("[" in raw_prefix and "]" in raw_prefix):
        raw_clusters = f"{raw_prefix}_clusters.csv"
        raw_info = f"{raw_prefix}_clustering_info.json"

        clusters_ok = _safe_write(raw_clusters, Path(args.out_clusters).read_text())
        info_ok = _safe_write(raw_info, info_text)

        if clusters_ok and info_ok:
            print("[INFO] Also wrote Nextflow-expected raw outputs:", flush=True)
            print(f"  - {raw_clusters}", flush=True)
            print(f"  - {raw_info}", flush=True)

    print(f"\n[SUCCESS] Saved:")
    print(f"  - {args.out_clusters}")
    print(f"  - {args.out_info}")
    print(f"{'='*60}\n", flush=True)


if __name__ == "__main__":
    main()