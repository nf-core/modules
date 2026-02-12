#!/usr/bin/env python3

import anndata as ad
import scvi
import yaml
from scipy.sparse import csr_matrix
from scvi.external import SCAR
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

scvi.settings.seed = 0


adata = ad.read_h5ad("${filtered}")
adata_unfiltered = ad.read_h5ad("${unfiltered}")

input_layer = "${input_layer ?: 'X'}"
output_layer = "${output_layer ?: 'scar'}"
max_epochs = "${max_epochs ?: ''}"
n_batch = int("${n_batch ?: 1}")

if input_layer != "X" and input_layer not in adata.layers:
    raise ValueError(f"Input layer {input_layer} not found in adata.layers")

SCAR.setup_anndata(adata, layer=None if input_layer == "X" else input_layer)
SCAR.get_ambient_profile(adata, adata_unfiltered, n_batch=n_batch)

vae = SCAR(adata)

# Prevent errors like https://discourse.scverse.org/t/scvi-21618-problem/2294
vae.train(
    max_epochs=int(max_epochs) if max_epochs else None,
    early_stopping=True,
    datasplitter_kwargs={"drop_last": True},
)

if output_layer == "X":
    adata.X = csr_matrix(vae.get_denoised_counts())
else:
    adata.layers[output_layer] = csr_matrix(vae.get_denoised_counts())

del adata.uns["_scvi_uuid"], adata.uns["_scvi_manager_uuid"]

adata.write_h5ad("${prefix}.h5ad")

# Versions
versions = {"${task.process}": {"scvi": scvi.__version__}}

with open("versions.yml", "w") as f:
    f.write(yaml.dump(versions))
