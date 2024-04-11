#!/usr/bin/env python

# This script is used to identify *.nf.test files for changed functions/processs/workflows/pipelines and *.nf-test files
# with changed dependencies, then return as a JSON list

import argparse
import json
import logging
import re
import yaml

from itertools import chain
from pathlib import Path
from git import Repo


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments and return an ArgumentParser object.

    Returns:
        argparse.ArgumentParser: The ArgumentParser object with the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Scan *.nf.test files for function/process/workflow name and return as a JSON list"
    )
    parser.add_argument(
        "-r",
        "--head_ref",
        required=True,
        help="Head reference branch (Source branch for a PR).",
    )
    parser.add_argument(
        "-b",
        "--base_ref",
        required=True,
        help="Base reference branch (Target branch for a PR).",
    )
    parser.add_argument(
        "-x",
        "--ignored_files",
        nargs="+",
        default=[
            ".git/*",
            ".gitpod.yml",
            ".prettierignore",
            ".prettierrc.yml",
            "*.md",
            "*.png",
            "modules.json",
            "pyproject.toml",
            "tower.yml",
        ],
        help="List of files or file substrings to ignore.",
    )
    parser.add_argument(
        "-i",
        "--include",
        type=Path,
        default=None,
        help="Path to an include file containing a YAML of key value pairs to include in changed files. I.e., return the current directory if an important file is changed.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level",
    )
    parser.add_argument(
        "-t",
        "--types",
        nargs="+",
        choices=["function", "process", "workflow", "pipeline"],
        default=["function", "process", "workflow", "pipeline"],
        help="Types of tests to include.",
    )
    return parser.parse_args()


def read_yaml_inverted(file_path: str) -> dict:
    """
    Read a YAML file and return its contents as a dictionary but reversed, i.e. the values become the keys and the keys become the values.

    Args:
        file_path (str): The path to the YAML file.

    Returns:
        dict: The contents of the YAML file as a dictionary inverted.
    """
    with open(file_path, "r") as f:
        data = yaml.safe_load(f)

    # Invert dictionary of lists into contents of lists are keys, values are the original keys
    # { "key": ["item1", "item2] } --> { "item1": "key", "item2": "key" }
    return {value: key for key, values in data.items() for value in values}


def find_changed_files(
    branch1: str,
    branch2: str,
    ignore: list[str],
) -> list[Path]:
    """
    Find all *.nf.tests that are associated with files that have been changed between two specified branches.

    Args:
        branch1 (str)      : The first branch being compared
        branch2 (str)      : The second branch being compared
        ignore  (list)     : List of files or file substrings to ignore.

    Returns:
        list: List of files matching the pattern *.nf.test that have changed between branch2 and branch1.
    """
    # create repo
    repo = Repo(".")
    # identify commit on branch1
    branch1_commit = repo.commit(branch1)
    # identify commit on branch2
    branch2_commit = repo.commit(branch2)
    # compare two branches
    diff_index = branch1_commit.diff(branch2_commit)

    # Start empty list of changed files
    changed_files = []

    # For every file that has changed between commits
    for file in diff_index:
        # Get pathlib.Path object
        filepath = Path(file.a_path)
        # If file does not match any in the ignore list, add containing directory to changed_files
        if not any(filepath.match(ignored_path) for ignored_path in ignore):
            changed_files.append(filepath)

    # Uniqueify the results before returning for efficiency
    return list(set(changed_files))


def detect_include_files(
    changed_files: list[Path], include_files: dict[str, str]
) -> list[Path]:
    """
    Detects the include files based on the changed files.

    Args:
        changed_files (list[Path]): List of paths to the changed files.
        include_files (dict[str, str]): Key-value pairs to return if a certain file has changed. If a file in a directory has changed, it points to a different directory.

    Returns:
        list[Path]: List of paths to representing the keys of the include_files dictionary, where a value matched a path in changed_files.
    """
    new_changed_files = []
    for filepath in changed_files:
        # If file is in the include_files, we return the key instead of the value
        for include_path, include_key in include_files.items():
            if filepath.match(include_path):
                new_changed_files.append(Path(include_key))
    return new_changed_files


def detect_files(paths: list[Path], suffix: str) -> list[Path]:
    """
    Detects and returns a list of nf-test files from the given list of changed files.

    Args:
        paths (list[Path]): A list of file paths to scan.
        suffix (str): File suffix to detect

    Returns:
        list[Path]: A list of nf-test file paths.
    """
    result: list[Path] = []

    for path in paths:
        # If Path is the exact nf-test file add to list:
        if path.match(suffix) and path.exists():
            result.append(path)
        # Else recursively search for nf-test files:
        elif path.is_dir():
            # Search the dir
            # e.g.
            # dir/
            # ├─ main.nf
            # ├─ main.nf.test
            for file in path.rglob(suffix):
                result.append(file)
        elif path.is_file():
            # Search the enclosing dir so files in the same dir can be found.
            # e.g.
            # dir/
            # ├─ main.nf
            # ├─ main.nf.test
            for file in path.parent.rglob(suffix):
                result.append(file)

    return result


def process_nf_test_files(files: list[Path]) -> list[str]:
    """
    Process the files and return lines that begin with 'workflow', 'process', or 'function' and have a single string afterwards.

    Args:
        files (list): List of files to process.

    Returns:
        list: List of lines that match the criteria.
    """
    result = []
    for file in files:
        with open(file, "r") as f:
            is_pipeline_test = True
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith(("workflow", "process", "function")):
                    words = line.split()
                    if len(words) == 2 and re.match(r'^".*"$', words[1]):
                        result.append(line)
                        is_pipeline_test = False

            # If no results included workflow, process or function
            # Add a dummy result to fill the 'pipeline' category
            if is_pipeline_test:
                result.append("pipeline 'PIPELINE'")

    return result


def convert_nf_test_files_to_test_types(
    lines: list[str], types: list[str] = ["function", "process", "workflow", "pipeline"]
) -> dict[str, list[str]]:
    """
    Generate a dictionary of function, process and workflow lists from the lines.

    Args:
        lines (list): List of lines to process.
        types (list): List of types to include.

    Returns:
        dict: Dictionary with function, process and workflow lists.
    """
    # Populate empty dict from types
    result: dict[str, list[str]] = {key: [] for key in types}

    for line in lines:
        words = line.split()
        if len(words) == 2 and re.match(r'^".*"$', words[1]):
            keyword = words[0]
            name = words[1].strip("'\"")  # Strip both single and double quotes
            if keyword in types:
                result[keyword].append(name)
        elif "run" in line:
            run_words = words[0].strip(")").split("(")
            if len(run_words) == 2 and re.match(r'^".*"$', run_words[1]):
                keyword = run_words[0]
                name = run_words[1].strip("'\"")  # Strip both single and double quotes
                if keyword in types:
                    result[keyword].append(name)
        elif "include {" in line and "include" in types:
            keyword = words[0]
            name = words[2].strip("'\"")  # Strip both single and double quotes
            if keyword in types:
                result[keyword].append(name)

    return result


def find_nf_tests_with_changed_dependencies(
    paths: list[Path], tags: list[str]
) -> list[Path]:
    """
    Find all *.nf.test files with changed dependencies
    (identified as modules loaded in the via setup { run("<tool>") } from a list of paths.

    Args:
        paths (list): List of directories or files to scan.
        tags (list): List of tags identified as having changes.

    Returns:
        list: List of *.nf.test files with changed dependencies.
    """

    result: list[Path] = []

    nf_test_files = detect_files(paths, "*.nf.test")

    # find nf-test files with changed dependencies
    for nf_test_file in nf_test_files:
        with open(nf_test_file, "r") as f:
            lines = f.readlines()
            # Get all tags from nf-test file
            # Make case insensitive with .casefold()
            tags_in_nf_test_file = [
                tag.casefold().replace("/", "_")
                for tag in convert_nf_test_files_to_test_types(lines, types=["run"])[
                    "run"
                ]
            ]
            # Check if tag in nf-test file appears in a tag.
            # Use .casefold() to be case insensitive
            if any(
                tag.casefold().replace("/", "_") in tags_in_nf_test_file for tag in tags
            ):
                result.append(nf_test_file)

    return result


def find_nf_files_with_changed_dependencies(
    paths: list[Path], tags: list[str]
) -> list[Path]:
    """
    Find all *.nf.test files with where the *.nf file uses changed dependencies
    (identified via include { <tool> }) in *.nf files from a list of paths.

    Args:
        paths (list): List of directories or files to scan.
        tags (list): List of tags identified as having changes.

    Returns:
        list: List of *.nf.test files from *.nf files with changed dependencies.
    """

    nf_files_w_changed_dependencies: list[Path] = []

    nf_files = detect_files(paths, "*.nf")

    # find nf files with changed dependencies
    for nf_file in nf_files:
        with open(nf_file, "r") as f:
            lines = f.readlines()
            # Get all include statements from nf file
            # Make case insensitive with .casefold()
            includes_in_nf_file = [
                tag.casefold().replace("/", "_")
                for tag in convert_nf_test_files_to_test_types(
                    lines, types=["include"]
                )["include"]
            ]
            # Check if include in nf file appears in a tag.
            # Use .casefold() to be case insensitive
            if any(
                tag.casefold().replace("/", "_") in includes_in_nf_file for tag in tags
            ):
                nf_files_w_changed_dependencies.append(nf_file)

    # find nf-test for nf files with changed dependencies
    nf_test_files_for_changed_dependencies = detect_files(
        nf_files_w_changed_dependencies, "*.nf.test"
    )

    return nf_test_files_for_changed_dependencies


if __name__ == "__main__":

    # Utility stuff
    args = parse_args()
    logging.basicConfig(level=args.log_level)

    # Parse nf-test files for target test tags
    changed_files = find_changed_files(args.head_ref, args.base_ref, args.ignored_files)

    # If an additional include YAML is added, we detect additional changed dirs to include
    if args.include:
        include_files = read_yaml_inverted(args.include)
        changed_files = changed_files + detect_include_files(
            changed_files, include_files
        )
    nf_test_files = detect_files(changed_files, "*.nf.test")
    lines = process_nf_test_files(nf_test_files)
    result = convert_nf_test_files_to_test_types(lines)

    # Get only relevant results (specified by -t)
    # Unique using a set
    target_results = list(
        {item for sublist in map(result.get, args.types) for item in sublist}
    )

    # Parse nf-test files to identify nf-tests containing "setup" with changed module/subworkflow/workflow
    nf_test_changed_setup = find_nf_tests_with_changed_dependencies(
        [Path("./modules/"), Path("./subworkflows/"), Path("./workflows/")],
        target_results,
    )

    # Parse *.nf files to identify nf-files containing include with changed module/subworkflow/workflow
    nf_files_changed_include = find_nf_files_with_changed_dependencies(
        [Path("./modules/"), Path("./subworkflows/"), Path("./workflows/")],
        target_results,
    )

    # Combine target nf-test files and nf-test files with changed dependencies
    all_nf_tests = [
        str(test_path)
        for test_path in set(
            nf_test_files + nf_test_changed_setup + nf_files_changed_include
        )
    ]

    # Print to stdout
    print(json.dumps(all_nf_tests))
