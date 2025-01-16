#!/usr/bin/env python3

import os
import platform

os.environ["MPLCONFIGDIR"] = "./tmp"

import anndata as ad
import scvi
import torch
from scvi.external import SOLO
from scvi.model import SCVI
from threadpoolctl import threadpool_limits

torch.set_float32_matmul_precision("medium")
scvi.settings.seed = 0

threadpool_limits(int("${task.cpus}"))
scvi.settings.num_threads = int("${task.cpus}")


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


def train_model(model):
    if "${task.ext.use_gpu}" == "true":
        model.to_device(0)

    model.train(max_epochs=int("${max_epochs}") if "${max_epochs}" else None, datasplitter_kwargs={"drop_last": True})


adata = ad.read_h5ad("${h5ad}")

if "${batch_key}":
    SCVI.setup_anndata(adata, batch_key="${batch_key}")
else:
    SCVI.setup_anndata(adata)

model = SCVI(adata)

train_model(model)

adata.obs["doublet"] = False

batches = adata.obs["batch"].unique() if "${batch_key}" else [0]
for batch in batches:
    model = SOLO.from_scvi_model(model, restrict_to_batch=batch if len(batches) > 1 else None)

    train_model(model)
    result = model.predict(False)

    doublets = result[result == "doublet"].index.tolist()
    adata.obs.loc[doublets, "doublet"] = True


df = adata.obs[["doublet"]]
df.columns = ["${prefix}"]
df.to_pickle("${prefix}.pkl")

adata = adata[~adata.obs["doublet"]].copy()
adata.obs.drop("doublet", axis=1, inplace=True)

del adata.uns["_scvi_manager_uuid"]
del adata.uns["_scvi_uuid"]

adata.write_h5ad("${prefix}.h5ad")

# Versions

versions = {
    "${task.process}": {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "scvi": scvi.__version__,
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
