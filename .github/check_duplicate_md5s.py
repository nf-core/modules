#!/usr/bin/env python

from rich import print
from rich.table import Table
import click
import glob
import os
import yaml


@click.command()
@click.option(
    "--min_dups",
    default=5,
    show_default=True,
    help="Minimum number of duplicates to report",
)
@click.option(
    "--search_dir",
    default=f"{os.path.dirname(__file__)}/../tests/**/test.yml",
    show_default=True,
    help="Glob directory pattern used to find test YAML files",
)
def find_duplicate_md5s(min_dups, search_dir):
    """
    Find duplicate file MD5 sums in test YAML files.
    """
    md5_filenames = {}
    md5_counts = {}

    # Loop through all files in tests/ called test.yml
    for test_yml in glob.glob(search_dir, recursive=True):
        # Open file and parse YAML
        with open(test_yml, "r") as fh:
            test_config = yaml.safe_load(fh)
            # Loop through tests and check for duplicate md5s
            for test in test_config:
                for test_file in test.get("files", []):
                    if "md5sum" in test_file:
                        md5 = test_file["md5sum"]
                        md5_filenames[md5] = md5_filenames.get(md5, []) + [
                            os.path.basename(test_file.get("path"))
                        ]
                        md5_counts[md5] = md5_counts.get(md5, 0) + 1

    # Set up rich table
    table = Table(title="Duplicate MD5s", row_styles=["dim", ""])
    table.add_column("MD5", style="cyan", no_wrap=True)
    table.add_column("Count", style="magenta", justify="right")
    table.add_column("Filenames", style="green")

    # Add rows - sort md5_counts by value
    for md5 in sorted(md5_counts, key=md5_counts.get):
        if md5_counts[md5] >= min_dups:
            table.add_row(md5, str(md5_counts[md5]), ", ".join(set(md5_filenames[md5])))

    print(table)


if __name__ == "__main__":
    find_duplicate_md5s()
