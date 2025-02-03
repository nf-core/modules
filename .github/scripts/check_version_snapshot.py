#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///

"""Check that versions are captured not hashed."""

import argparse
import re
from pathlib import Path
from typing import Optional, Sequence


def check_snapshotting(lines: list[str], path: Path) -> None:
    """Check each line for proper snapshotting."""
    good_pattern = re.compile(
        r'snapshot\(path\(process\.out\.versions\.get\(0\)\)\.yaml\)\.match\("versions"\)|path\(process\.out\.versions\.get\(0\)\)\.yaml'
    )
    bad_pattern = re.compile(r"process\.out\.versions|snapshot\(process\.out\)\.match\(\)")

    for line_number, line in enumerate(lines, start=1):
        if bad_pattern.search(line) and not good_pattern.search(line):
            print(f"Improper snapshotting in {path} at line {line_number}: {line.strip()}")


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Check that versions are captured not hashed."""
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    args = parser.parse_args(argv)
    for path in args.paths:
        with path.open() as f:
            lines = f.readlines()
            check_snapshotting(lines, path)


if __name__ == "__main__":
    main()
