#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "ruamel.yaml",
# ]
# ///

"""Sort dependencies in conda environment files."""
# Test with
# uv run --with pytest --with "ruamel.yaml" -- pytest -v .github/scripts/conda_env_sorter.py

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

import ruamel.yaml

# Add pytest imports conditionally
if "pytest" in sys.modules:
    import pytest

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

        # Sort channels if they exist
        if "channels" in doc:
            doc["channels"].sort(key=str)

        # Sort dependencies
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

# Pytest tests (only loaded when running pytest)
if "pytest" in sys.modules:

    @pytest.mark.parametrize(
        "input_content,expected",
        [
            # Test basic sorting
            ("dependencies:\n  - zlib\n  - python\n", ["python", "zlib"]),
            # Test dict sorting
            ("dependencies:\n  - pip:\n    - b\n    - a\n  - python\n", ["python", {"pip": ["a", "b"]}]),
            # Test existing headers
            ("---\n# yaml-language-server: $schema=...\ndependencies:\n  - b\n  - a\n", ["a", "b"]),
        ],
    )
    def test_conda_sorter(tmp_path, input_content, expected):
        test_file = tmp_path / "environment.yml"
        test_file.write_text(input_content)

        # Run our sorter on the test file
        main([str(test_file)])

        # Read back the sorted file
        result = test_file.read_text()

        # Check schema headers are present
        assert result.startswith("---\n# yaml-language-server: $schema=")

        # Parse the sorted content (skip first 2 header lines)
        parsed = yaml.load("".join(result.splitlines(True)[2:]))

        # Compare the actual dependencies structure
        assert parsed["dependencies"] == expected

    def test_invalid_file(tmp_path):
        test_file = tmp_path / "bad.yml"
        test_file.write_text("invalid: yaml: here")

        with pytest.raises(ruamel.yaml.scanner.ScannerError):
            main([str(test_file)])

    def test_full_conda_env_with_channels(tmp_path):
        """Test sorting a full conda environment file with channels and namespaced dependencies."""
        test_file = tmp_path / "full_env.yml"
        test_file.write_text(input_content.strip())

        # Run our sorter on the test file
        main([str(test_file)])

        # Read back the sorted file
        result = test_file.read_text()

        # Check that result matches the expected output
        assert result.strip() == output_content.strip()

# Test input and output for full conda environment with channels
input_content = """
---
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
channels:
  - conda-forge
  - bioconda
dependencies:
  - bioconda::ngscheckmate=1.0.1
  - bioconda::bcftools=1.21
"""

# Should be sorted because it should sort the channels first, then the dependency names
output_content = """
---
# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
channels:
  - bioconda
  - conda-forge
dependencies:
  - bioconda::bcftools=1.21
  - bioconda::ngscheckmate=1.0.1
"""
