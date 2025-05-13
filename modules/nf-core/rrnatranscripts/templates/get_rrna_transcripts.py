#!/usr/bin/env python3

import logging
import platform
import sys
from pathlib import Path

# Configure logging
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_rrna_intervals(gtf: str, rrna_transcripts: str):
    """
    Get lines containing ``#`` or ``gene_type rRNA`` or ```` or ``gene_type rRNA_pseudogene`` or ``gene_type MT_rRNA``
    Create output file

    Args:
        file_in (pathlib.Path): The given GTF file.
        file_out (pathlib.Path): Where the ribosomal RNA GTF file should
            be created; always in GTF format.
    """
    patterns = {
        "#",
        'transcript_biotype "Mt_rRNA"',
        'transcript_biotype "rRNA"',
        'transcript_biotype "rRNA_pseudogene"',
    }
    line_starts = {"MT", "1", "2", "3", "4", "5", "6", "7", "8", "9"}
    out_lines = []
    path_gtf = Path(gtf)
    path_rrna_transcripts = Path(rrna_transcripts)
    if not path_gtf.is_file():
        logger.error(f"The given input file {gtf} was not found!")
        sys.exit(2)
    with path_gtf.open() as f:
        data = f.readlines()
        for line in data:
            for pattern in patterns:
                if pattern in line:
                    for line_start in line_starts:
                        if line.startswith(line_start):
                            out_lines.append(line)
    if out_lines != []:
        with path_rrna_transcripts.open(mode="w") as out_file:
            out_file.writelines(out_lines)


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


if __name__ == "__main__":
    if "${task.ext.prefix}" != "null":
        prefix = "${task.ext.prefix}."
    else:
        prefix = "${task.ext.gtf}."

    if not get_rrna_intervals("$gtf", f"{prefix}_rrna_intervals.gtf"):
        logging.error("Failed to extract rrna transcipts.")

    # Write the versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))
