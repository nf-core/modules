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
        if doc is None:
            raise ruamel.yaml.scanner.ScannerError("Empty YAML file")

        # Sort channels if they exist
        if "channels" in doc:
            # Define channel priority order
            channel_order = {
                "conda-forge": 0,
                "bioconda": 1,
            }

            # Sort channels based on priority order, then alphabetically within same priority
            doc["channels"].sort(key=lambda x: (channel_order.get(x, 2), str(x)))

        # Sort dependencies if they exist
        if "dependencies" in doc:
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
            # Test channel sorting
            (
                "channels:\n  - conda-forge\n  - bioconda\ndependencies:\n  - python\n",
                {"channels": ["conda-forge", "bioconda"], "dependencies": ["python"]},
            ),
            # Test channel sorting with additional channels
            (
                "channels:\n  - bioconda\n  - conda-forge\n  - defaults\n  - r\n",
                {"channels": ["conda-forge", "bioconda", "defaults", "r"]},
            ),
            # Test namespaced dependencies
            (
                "dependencies:\n  - bioconda::ngscheckmate=1.0.1\n  - bioconda::bcftools=1.21\n",
                ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1"],
            ),
            # Test mixed dependencies
            (
                "dependencies:\n  - bioconda::ngscheckmate=1.0.1\n  - python\n  - bioconda::bcftools=1.21\n",
                ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1", "python"],
            ),
            # Test full environment with channels and namespaced dependencies
            (
                """
                channels:
                    - conda-forge
                    - bioconda
                dependencies:
                    - bioconda::ngscheckmate=1.0.1
                    - bioconda::bcftools=1.21
                """,
                {
                    "channels": ["conda-forge", "bioconda"],
                    "dependencies": ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1"],
                },
            ),
        ],
        ids=[
            "basic_dependency_sorting",
            "dict_dependency_sorting",
            "existing_headers",
            "channel_sorting",
            "channel_sorting_with_additional_channels",
            "namespaced_dependencies",
            "mixed_dependencies",
            "full_environment",
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
        if isinstance(expected, list):
            assert parsed["dependencies"] == expected
        else:
            # For comparing dictionaries, only compare the keys that are in the expected dictionary
            for key, value in expected.items():
                assert key in parsed
                assert parsed[key] == value

    def test_invalid_file(tmp_path):
        test_file = tmp_path / "bad.yml"
        test_file.write_text("invalid: yaml: here")

        with pytest.raises(ruamel.yaml.scanner.ScannerError):
            main([str(test_file)])

    def test_empty_file(tmp_path):
        """Test handling of empty files."""
        test_file = tmp_path / "empty.yml"
        test_file.write_text("")

        with pytest.raises(ruamel.yaml.scanner.ScannerError):
            main([str(test_file)])

    def test_missing_dependencies(tmp_path):
        """Test handling of files without dependencies section."""
        test_file = tmp_path / "no_deps.yml"
        test_file.write_text("channels:\n  - conda-forge\n")

        # Run without error now that we handle missing dependencies
        main([str(test_file)])

        # Read back and verify channel is preserved
        result = test_file.read_text()
        parsed = yaml.load("".join(result.splitlines(True)[2:]))
        assert "channels" in parsed
        assert parsed["channels"] == ["conda-forge"]
        assert "dependencies" not in parsed
