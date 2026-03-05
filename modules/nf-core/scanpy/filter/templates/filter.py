#!/usr/bin/env python3

import os

# These are needed to prevent errors during import of scanpy
# when using singularity/apptainer
os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import platform

import scanpy as sc
import yaml
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

adata = sc.read_h5ad("${h5ad}")
prefix = "${prefix}"
symbol_col = "${symbol_col}"

if symbol_col != "index" and symbol_col not in adata.var.columns:
    raise ValueError(f"Column '{symbol_col}' not found in adata.var: '{adata.var.columns}'")

symbol_series = adata.var.index if symbol_col == "index" else adata.var[symbol_col]
adata.var["mt"] = symbol_series.str.lower().str.startswith("mt-")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < int("${max_mito_percentage}"), :].copy()

sc.pp.filter_cells(adata, min_counts=int("${min_counts_cell}"))
sc.pp.filter_genes(adata, min_counts=int("${min_counts_gene}"))

sc.pp.filter_cells(adata, min_genes=int("${min_genes}"))
sc.pp.filter_genes(adata, min_cells=int("${min_cells}"))

adata.write_h5ad(f"{prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "scanpy": sc.__version__}}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
