#!/usr/bin/env python

import yaml
from metaspace import SMInstance

with open("${config}") as f:
    config = yaml.safe_load(f)

sm = SMInstance(api_key=config['User_info']['API_key'])

ds_name = config['Dataset_name']

Annot_settings = config['Annotation_settings']

metadata = {k: v for k, v in config.items() if k not in ['User_info','Dataset_name','Annotation_settings']}
metadata = {'Data_Type': 'Imaging MS', **metadata}

dbs = [tuple(item) for item in Annot_settings["Metabolite_database"]]

if len(Annot_settings['Adducts']) == 0:
    adduct_list = None
else:
    adduct_list = Annot_settings['Adducts']

if len(Annot_settings['Chemical_modifications']) == 0:
    chem_modif = None
else:
    chem_modif = Annot_settings['Chemical_modifications']

scoring_model_id_map = {
    'MSM': 2,
    'METASPACE_ML_Animal': 3,
    'METASPACE_ML_Plant': 4
}

scoring_model_id = scoring_model_id_map[Annot_settings['Analysis_version']]

# Submit dataset

ds_id = sm.submit_dataset(
    imzml_fn = "${imzml}",
    ibd_fn = "${ibd}",
    name = ds_name,
    metadata = metadata,
    is_public = False,
    databases = dbs,
    adducts = adduct_list,
    chem_mods = chem_modif,
    scoring_model = scoring_model_id,
    ppm = float(Annot_settings['ppm'])
)

with open("ds_id.txt", 'w') as f:
    f.write(ds_id)

# Versions
versions = {"${task.process}": {"python": "3.11",
                                "metaspace": "2.0.9",
                                "yaml": "6.0.2"}}

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
