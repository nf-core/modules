#!/usr/bin/env python3

"""
Author:
    Leon Hafner

Copyright (c) 2024, Leon Hafner. All rights reserved.

License: MIT License
"""

import platform

import anndata as ad


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


adata = ad.read_h5ad("${h5ad}", backed="r")

size_functions = {
    "cells": lambda: adata.n_obs,
    "genes": lambda: adata.n_vars,
    "obs": lambda: adata.n_obs,
    "var": lambda: adata.n_vars,
}

size_type = "${size_type}".lower()
if size_type not in size_functions:
    raise ValueError(f'Size type must be one of {', '.join(size_functions.keys())}.')

size = size_functions[size_type]()

with open("${prefix}.txt", "w") as f:
    f.write(str(size))


# Versions
versions = {"${task.process}": {"python": platform.python_version(), "anndata": ad.__version__}}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))
