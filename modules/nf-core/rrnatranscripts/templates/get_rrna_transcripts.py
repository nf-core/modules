#!/usr/bin/env python3

import argparse
import logging
import platform
import sys
from pathlib import Path


def get_rrna_intervals(file_in, file_out):
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
    with file_in.open() as f:
        data = f.readlines()
        for line in data:
            for pattern in patterns:
                if pattern in line:
                    for line_start in line_starts:
                        if line.startswith(line_start):
                            out_lines.append(line)

    with file_out.open(mode="w") as out_file:
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


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract ribosomal RNA intervals from a gtf file.",
        epilog="Example: python get_rrna_transcripts.py <filename.gtf> <intervals.gtf>",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input in GTF format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output intervals in GTF format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    get_rrna_intervals(args.file_in, args.file_out)

    # Write the versions
    versions_this_module = {}
    versions_this_module["${task.process}"] = {"python": platform.python_version()}
    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions_this_module))


if __name__ == "__main__":
    sys.exit(main())
