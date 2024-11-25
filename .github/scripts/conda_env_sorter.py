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


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Sort dependencies in conda environment files."""
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)
    for path in args.paths:
        doc = yaml.load(path)
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

        yaml.dump(doc, path)


if __name__ == "__main__":
    main()
