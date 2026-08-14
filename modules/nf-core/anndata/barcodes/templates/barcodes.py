#!/usr/bin/env python3
"""
Filter an AnnData object based on a list of barcodes.

This script reads a CSV file containing cell barcodes and an HDF5 file with an AnnData object,
then subsets the AnnData object to include only cells whose barcodes are present in the provided list.
The filtered AnnData object is written to a new HDF5 file.

Inputs:
    - ${barcodes}: CSV file with a single column of cell barcodes (one barcode per line, no header).
    - ${h5ad}: HDF5 file containing an AnnData object.

Outputs:
    - ${prefix}.h5ad: HDF5 file containing the filtered AnnData object.
    - versions.yml: YAML file with tool version information.
"""
import platform

import anndata as ad
import pandas as pd
import yaml

df = pd.read_csv("${barcodes}", header=None)
adata = ad.read_h5ad("${h5ad}")

adata = adata[df[0].values]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {"python": platform.python_version(), "anndata": ad.__version__, "pandas": pd.__version__}
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
