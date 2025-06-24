#!/usr/bin/env python3

import platform

import anndata as ad
import cellbender
from cellbender.remove_background.downstream import load_anndata_from_input_and_output

adata = ad.read_h5ad("${filtered}")

adata_cellbender = load_anndata_from_input_and_output("${unfiltered}", "${cellbender_h5}", analyzed_barcodes_only=False)

adata_cellbender = adata_cellbender[adata.obs_names]

if "${output_layer}" == "X":
    adata.X = adata_cellbender.layers["cellbender"]
else:
    adata.layers["${output_layer}"] = adata_cellbender.layers["cellbender"]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {"${task.process}": {"python": platform.python_version(), "cellbender": cellbender.__version__}}
