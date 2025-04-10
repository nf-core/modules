#!/usr/bin/env python3

import platform

import anndata as ad
import scvi
from scipy.sparse import csr_matrix
from scvi.external import SCAR
from threadpoolctl import threadpool_limits

threadpool_limits(int("${task.cpus}"))

scvi.settings.seed = 0


def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "    " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str


adata = ad.read_h5ad("${filtered}")
adata_unfiltered = ad.read_h5ad("${unfiltered}")

SCAR.setup_anndata(adata, layer=None if "${input_layer}" == "X" else "${input_layer}")
SCAR.get_ambient_profile(adata, adata_unfiltered)

vae = SCAR(adata)

# Prevent errors like https://discourse.scverse.org/t/scvi-21618-problem/2294
vae.train(
    max_epochs=int("${max_epochs}") if "${max_epochs}" else None,
    early_stopping=True,
    datasplitter_kwargs={"drop_last": True},
)

if "${output_layer}" == "X":
    adata.X = csr_matrix(vae.get_denoised_counts())
else:
    adata.layers["${output_layer}"] = csr_matrix(vae.get_denoised_counts())

del adata.uns["_scvi_uuid"], adata.uns["_scvi_manager_uuid"]

adata.write_h5ad("${prefix}.h5ad")

# Versions
versions = {"${task.process}": {"python": platform.python_version(), "scvi": scvi.__version__}}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
