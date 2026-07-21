#!/usr/bin/env python3

import json
import platform
import sklearn
import yaml
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, DBSCAN

PC_COL_RE = re.compile('[Pp][Cc][0-9]+', re.IGNORECASE)


def convert_eigenvec_to_tsv(eigenvec_path, out_pca, id_mode='iid'):
    rows = []
    n_pcs = 0
    mode = None

    with eigenvec_path.open('r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if parts[0].startswith('#'):
                header = [p.lstrip('#') for p in parts]
                if len(header) >= 2 and header[0].upper() == 'FID' and header[1].upper() == 'IID':
                    mode = 'fid_iid'
                elif header[0].upper() == 'IID':
                    mode = 'iid_only'
                continue
            if mode is None:
                try:
                    float(parts[1])
                    mode = 'iid_only'
                except (ValueError, IndexError):
                    mode = 'fid_iid'
            if mode == 'fid_iid':
                if len(parts) < 3:
                    continue
                fid = parts[0]
                iid = parts[1]
                pcs = parts[2:]
                sample_id = iid if id_mode == 'iid' else f'{fid}:{iid}'
            elif mode == 'iid_only':
                if len(parts) < 2:
                    continue
                iid = parts[0]
                pcs = parts[1:]
                sample_id = iid
            else:
                raise ValueError(f'Unrecognized eigenvec format in {eigenvec_path}')
            if n_pcs == 0:
                n_pcs = len(pcs)
            rows.append((sample_id, pcs))

    if not rows:
        raise ValueError(f'No valid data found in {eigenvec_path}')

    header = ['sample_id'] + [f'PC{i+1}' for i in range(n_pcs)]
    with out_pca.open('w') as fh:
        fh.write('\\t'.join(header) + '\\n')
        for sample_id, pcs in rows:
            fh.write(sample_id + '\\t' + '\\t'.join(pcs) + '\\n')

    print(f'[INFO] Converted {len(rows)} samples with {n_pcs} PCs -> {out_pca}')
    return n_pcs


def read_table_robust(path):
    df = pd.read_csv(path, sep='\\t', dtype=str)
    print(f'[DEBUG] Initial read: {df.shape[0]} rows x {df.shape[1]} cols', flush=True)
    col_names_upper = set(str(c).upper() for c in df.columns)

    def is_header_row(row):
        row_values_upper = [str(v).upper() for v in row.values]
        overlap = sum(1 for v in row_values_upper if v in col_names_upper)
        if overlap >= 3:
            return True
        header_keywords = {'FID', 'IID', 'PC1', 'PC2', 'PC3'}
        if sum(1 for v in row_values_upper if v in header_keywords) >= 2:
            return True
        return False

    bad_rows = df.apply(is_header_row, axis=1)
    if bad_rows.any():
        n_bad = int(bad_rows.sum())
        print(f'[INFO] Removed {n_bad} duplicate header row(s)', flush=True)
        df = df[~bad_rows].copy().reset_index(drop=True)

    print(f'[INFO] After cleanup: {df.shape[0]} rows x {df.shape[1]} cols', flush=True)
    return df


def build_sample_id(df):
    cols = list(df.columns)
    if 'sample_id' in df.columns:
        return df['sample_id'].astype(str), df.drop(columns=['sample_id'])
    iid_candidates = [c for c in cols if str(c).upper() == 'IID']
    if iid_candidates:
        iid = iid_candidates[0]
        return df[iid].astype(str), df.drop(columns=[iid])
    fid_candidates = [c for c in cols if str(c).upper() == 'FID']
    if fid_candidates and iid_candidates:
        fid = fid_candidates[0]
        iid = iid_candidates[0]
        sample_ids = df[iid].astype(str)
        return sample_ids, df.drop(columns=[c for c in [fid, iid] if c in df.columns])
    pc_cols = [c for c in cols if PC_COL_RE.match(str(c))]
    non_pc_cols = [c for c in cols if c not in pc_cols]
    if non_pc_cols:
        id_col = non_pc_cols[0]
        return df[id_col].astype(str), df.drop(columns=[id_col])
    return pd.Series([f'sample_{i}' for i in range(len(df))], index=df.index), df


def main():
    prefix = '${meta.id}'

    pca_tsv = Path(f'{prefix}_pca_scores.tsv')
    convert_eigenvec_to_tsv(Path('${eigenvec}'), pca_tsv, 'iid')

    df = read_table_robust(str(pca_tsv))
    sample_ids, df_feats = build_sample_id(df)

    pc_cols = [c for c in df_feats.columns if PC_COL_RE.match(str(c))]
    if not pc_cols:
        raise ValueError('No PC columns found in input')

    X = df_feats[pc_cols].apply(pd.to_numeric, errors='coerce').values
    if np.isnan(X).any():
        raise ValueError('NaN values detected in PCA data')

    print(f'[INFO] Loaded {X.shape[0]} samples x {X.shape[1]} principal components', flush=True)

    if '${algorithm}' == 'kmeans':
        model = KMeans(n_clusters=${n_clusters}, init='random', n_init=100, random_state=42)
        labels = model.fit_predict(X)
        info = {'algorithm': 'kmeans', 'k': ${n_clusters}, 'inertia': float(model.inertia_)}
    else:
        model = DBSCAN(eps=${dbscan_eps}, min_samples=${dbscan_min_samples})
        labels = model.fit_predict(X)
        n_found = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = int(np.sum(labels == -1))
        info = {
            'algorithm': 'dbscan',
            'eps': ${dbscan_eps},
            'min_samples': ${dbscan_min_samples},
            'n_clusters_found': int(n_found),
            'n_noise': n_noise
        }

    out_clusters = f'{prefix}_clusters.csv'
    out_info = f'{prefix}_clustering_info.json'

    pd.DataFrame({'sample_id': sample_ids.astype(str), 'cluster': labels}).to_csv(out_clusters, index=False)
    info.update({
        'n_samples': int(X.shape[0]),
        'n_features': int(X.shape[1]),
        'feature_names': pc_cols,
        'input_file': Path('${eigenvec}').name
    })
    Path(out_info).write_text(json.dumps(info, indent=2))

    print('[SUCCESS] Clustering completed:')
    print(f'   -> Clusters : {out_clusters}')
    print(f'   -> Info     : {out_info}')


    versions = {
        'CUSTOM_CLUSTERING': {
            'python': platform.python_version(),
            'scikit-learn': sklearn.__version__,
            'pandas': pd.__version__,
            'numpy': np.__version__,
        }
    }
    with open('versions.yml', 'w') as fh:
        fh.write(yaml.dump(versions, default_flow_style=False))

main()
