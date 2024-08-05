#!/usr/bin/env python3

import scvi
import anndata as ad
from scvi.external import SCAR
import platform

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))


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


adata = ad.read_h5ad("${h5ad}")
adata_raw = ad.read_h5ad("${raw}")

# TODO: Find out why the batch_key='batch' argument causes an error.
SCAR.setup_anndata(adata)
SCAR.get_ambient_profile(adata, adata_raw)

vae = SCAR(adata)

# The training fails if the number of entries in a minibatch is 1.
# This happens if the total number of entries modulo the batch size is 1.
# We therefore decrease the batch size until this is not the case.
batch_size = 128
worked = False
while not worked:
    try:
        vae.train(
            batch_size=batch_size,
            early_stopping=True
        )
        worked = True
    except ValueError as e:
        batch_size -= 1

        if batch_size < 125:
            raise e

adata.layers["ambient"] = vae.get_denoised_counts()

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scvi": scvi.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
