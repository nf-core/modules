#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "ruamel.yaml",
# ]
# ///

"""Sort dependencies in conda environment files."""

import argparse
from pathlib import Path
from typing import Optional, Sequence

import ruamel.yaml

yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=2, offset=2)  # Set indentation to 2 spaces


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Sort dependencies in conda environment files."""
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)
    for path in args.paths:
        # Read the entire file content
        with path.open() as f:
            lines = f.readlines()

        # Define the schema lines to be added if missing
        schema_lines = [
            "---\n",
            "# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json\n",
        ]

        # Check if the first two lines match the expected schema lines
        if lines[:2] == schema_lines:
            content = "".join(lines[2:])  # Skip schema lines when reading content
        else:
            content = "".join(lines)  # Use all content if no schema lines present

        # Parse the YAML content
        doc = yaml.load(content)
        dicts = []
        others = []

        for term in doc["dependencies"]:
            if isinstance(term, dict):
                dicts.append(term)
            else:
                others.append(term)
        others.sort(key=str)
        for dict_term in dicts:
            for value in dict_term.values():
                if isinstance(value, list):
                    value.sort(key=str)
        dicts.sort(key=str)
        doc["dependencies"].clear()
        doc["dependencies"].extend(others)
        doc["dependencies"].extend(dicts)

        # Write back to file with headers
        with path.open("w") as f:
            # Always write schema lines first
            f.writelines(schema_lines)
            # Then dump the sorted YAML
            yaml.dump(doc, f)


if __name__ == "__main__":
    main()
