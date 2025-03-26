#!/usr/bin/env python

import os
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'

import metaspace_converter as mc
from metaspace import SMInstance
import sys

def str_to_bool(value):
    if value.lower() in ['true', '1', 't', 'y', 'yes']:
        return True
    elif value.lower() in ['false', '0', 'f', 'n', 'no']:
        return False
    else:
        raise ValueError('Boolean value expected.')

sm = SMInstance()
ds = sm.dataset(id = "${ds_id}")

annot_dbs = [(db.name,db.version) for db in ds.database_details]
annot_lengths = [len(ds.annotations(database = i,
                                    fdr = float("${fdr}"))) for i in annot_dbs]
max_annot_idx = annot_lengths.index(max(annot_lengths))
input_db = ("${database_name}", "${database_version}")
if input_db not in annot_dbs:
    db = annot_dbs[max_annot_idx]
else:
    db = input_db

adata = mc.metaspace_to_anndata(
    dataset=ds,
    use_tic = str_to_bool("${use_tic}"),
    metadata_as_obs = str_to_bool("${metadata_as_obs}"),
    fdr=float("${fdr}"),
    database=db
)

sdata = mc.metaspace_to_spatialdata(
    dataset=ds,
    use_tic = str_to_bool("${use_tic}"),
    metadata_as_obs = str_to_bool("${metadata_as_obs}"),
    fdr=float("${fdr}"),
    database=db)

adata.write_h5ad("AnnData_${ds_id}.h5ad")
sdata.write("SpatialData_${ds_id}.zarr", overwrite=True)

# Versions
versions = {"${task.process}": {"python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
                                "metaspace_converter": "1.1.1"}}

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

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
