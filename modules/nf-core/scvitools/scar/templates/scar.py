#!/usr/bin/env python3

import scvi
import anndata as ad
from scvi.external import SCAR
import platform

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
        spaces = "  " * indent
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

# The training fails if the number of entries in a minibatch is 1.
# This happens if the total number of cells modulo the batch size is 1.
batch_size = [size for size in range(128, 120, -1) if adata.n_obs % size != 1][0]
vae.train(batch_size=batch_size, early_stopping=True)

if "${output_layer}" == "X":
    adata.X = vae.get_denoised_counts()
else:
    adata.layers["${output_layer}"] = vae.get_denoised_counts()

del adata.uns["_scvi_uuid"], adata.uns["_scvi_manager_uuid"]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "scvi": scvi.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
