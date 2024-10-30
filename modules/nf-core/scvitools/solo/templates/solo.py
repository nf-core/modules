#!/usr/bin/env python3

import platform
import random

import numpy as np
import anndata as ad
import scvi
import torch
from scvi.external import SOLO
from scvi.model import SCVI
from threadpoolctl import threadpool_limits

torch.set_float32_matmul_precision("medium")
scvi.settings.seed = 0
torch.manual_seed(0)
np.random.seed(0)
random.seed(0)
torch.use_deterministic_algorithms(True)

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
    def generate_batch_sizes():
        attempts = 0
        while True:
            yield 128 + 32 * attempts
            attempts += 1

    if "${task.ext.use_gpu}" == "true":
        model.to_device(0)

    for batch_size in generate_batch_sizes():
        try:
            model.train(batch_size=batch_size, max_epochs=int("${max_epochs}"))
            break
        except Exception as e:
            print(f"Failed with batch size {batch_size}: {e}")


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
