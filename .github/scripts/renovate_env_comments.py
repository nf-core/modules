#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Optional, Sequence

# Define a regex pattern to match dependencies
dependency_pattern = re.compile(r"^\s*-\s*(?P<channel>[\w-]+)::(?P<package>[\w-]+)=?(?P<version>[\d\.]*)")


# Function to add or verify renovate comments
def process_environment_file(file_path):
    with open(file_path) as file:
        lines = file.readlines()

    updated_lines = []
    for line in lines:
        match = dependency_pattern.match(line)
        if match:
            channel = match.group("channel")
            package = match.group("package")
            version = match.group("version")
            comment = f"  # renovate: datasource=conda depName={channel}/{package}\n"

            # Check if the comment is already present
            if comment not in updated_lines:
                updated_lines.append(comment)

        updated_lines.append(line)

    # Write the updated lines back to the file
    with open(file_path, "w") as file:
        file.writelines(updated_lines)


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Add Renovate comments to dependencies in conda environment files."""
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)
    for path in args.paths:
        print(path)
        process_environment_file(path)


if __name__ == "__main__":
    main()
