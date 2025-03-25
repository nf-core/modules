#!/usr/bin/env python

import os
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'

import metaspace_converter as mc
from metaspace import SMInstance
import argparse

def str_to_bool(value):
    if value.lower() in ['true', '1', 't', 'y', 'yes']:
        return True
    elif value.lower() in ['false', '0', 'f', 'n', 'no']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# parser = argparse.ArgumentParser(description="Convert a METASPACE dataset to AnnData and SpatialData")

# parser.add_argument('--fdr', type=int, help='FDR threshold')
# parser.add_argument('--use_tic', type=str_to_bool, help='Normalize by TIC')
# parser.add_argument('--metadata_as_obs', type=str_to_bool, help='Add metaspace metadata as obs')
# parser.add_argument('--database_name', type=str, help='Database name')
# parser.add_argument('--database_version', type=str, help='Database version')

# # Parse the arguments
# args = parser.parse_args()

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
versions = {"${task.process}": {"python": "3.11", "metaspace_converter": "1.1.1"}}

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
