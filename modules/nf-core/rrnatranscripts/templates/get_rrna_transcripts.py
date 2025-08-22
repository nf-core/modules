#!/usr/bin/env python3

"""
Author:
    Annick Renevey

Copyright (c) 2025, Annick Renevey. All rights reserved.

License: GPL-3 License
"""

import platform
from pathlib import Path

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

def get_rrna_intervals(file_in: Path, file_out: Path):
    """
    Extracts rRNA lines from a GTF and writes them out.
    """
    patterns = {
        "#",
        'transcript_biotype "Mt_rRNA"',
        'transcript_biotype "rRNA"',
        'transcript_biotype "rRNA_pseudogene"',
    }
    line_starts = {"MT", "1", "2", "3", "4", "5", "6", "7", "8", "9"}
    out_lines = []

    with file_in.open() as f:
        for line in f:
            if any(p in line for p in patterns) and any(line.startswith(ls) for ls in line_starts):
                out_lines.append(line)

    file_out.parent.mkdir(parents=True, exist_ok=True)
    with file_out.open("w") as out:
        out.writelines(out_lines)

# Main
gtf_path  = Path("${gtf}")
prefix    = "${prefix}"
out_gtf   = Path(f"{prefix}_rrna_intervals.gtf")

get_rrna_intervals(gtf_path, out_gtf)

# Versions
versions = {
    "${task.process}": {
        "python": platform.python_version()
    }
}

with open("versions.yml", "w") as v:
    v.write(format_yaml_like(versions))
