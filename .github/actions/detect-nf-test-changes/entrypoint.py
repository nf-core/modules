#!/usr/bin/env python

# This script is used to identify *.nf.test files for changed functions/processs/workflows/pipelines and *.nf-test files
# with changed dependencies, then return as a JSON list

import argparse
import fnmatch
import json
import logging
import os
from enum import Enum
from pathlib import Path

import yaml
from git import Repo


class TestTargetType(Enum):
    """
    Represents the type of test target.
    """

    FUNCTION: str = "function"
    PROCESS: str = "process"
    WORKFLOW: str = "workflow"
    PIPELINE: str = "pipeline"


class NextflowFile:
    """
    Represents a Nextflow file.

    Attributes:
        path (str): The path to the Nextflow file.
        lines (list[str]): The lines of the Nextflow file.
        includes (list[str]): A list of imported Nextflow workflows, processes or functions
    """

    def __init__(self, path):
        self.path = path
        logging.debug(f"Reading {self.path}")
        self.lines = self.read_nf_file()
        self.includes = self.find_include_statements()

    def read_nf_file(self) -> list[str]:
        """
        Read the Nextflow file and return its lines as a list.

        Returns:
            list[str]: The lines of the Nextflow file.
        """
        with open(self.path, "r") as f:
            return f.readlines()

    def find_include_statements(self) -> list[str]:
        """
        Find all include statements in the Nextflow file.

        Returns:
            list[str]: The include statements found in the Nextflow file.
        """
        result = []
        for line in self.lines:
            if "include {" in line:
                # Handle multi-line include statements
                if "{" in line and "}" not in line:
                    # Keep reading lines until we find the closing brace
                    full_include = line
                    line_idx = self.lines.index(line) + 1
                    while (
                        line_idx < len(self.lines) and "}" not in self.lines[line_idx]
                    ):
                        full_include += self.lines[line_idx]
                        line_idx += 1
                    if line_idx < len(self.lines):
                        full_include += self.lines[line_idx]
                    line = full_include

                # Extract the 'from' part which contains the path
                dependency = line.split()[2].strip("'\"").replace("/", "_").casefold()
                result.append(dependency)
        return result


class NfTest:
    """
    Represents a nf-test file.

    Attributes:
        test_path (Path): The path to the test file.
        test_dir (Path): The directory containing the test file.
        lines (list[str]): The lines of the test file.
        test_type (TestTargetType): The type of the test.
        test_name (str): The name of the test.
        run_statements (list[str]): The run statements in the test.
        config_file (Path | None): The path to the configuration file, or None if it doesn't exist.
        nextflow_path (Path): The path to the Nextflow script.
        nextflow (NextflowFile): The NextflowFile object representing the Nextflow script.
        root_path (Path): The common path between the Nextflow script and the test file.
        dependencies (list[str]): The dependencies of the test.
        tags (list[str]): The tags of the test.
    """

    def __init__(self, path, repo: Path = Path(".")):
        self.test_path = path
        self.repo = repo
        self.populate_attributes()

    def __str__(self):
        return f"Test: {self.test_name}, Type: {self.test_type}, Path: {self.test_path}"

    def read_test_file(self, readpath: Path) -> list[str]:
        with open(readpath, "r") as f:
            return f.readlines()

    def populate_attributes(self):
        self.test_dir = self.test_path.parent
        self.lines = self.read_test_file(self.test_path)
        self.test_type, self.test_name = self.find_test_type()
        self.run_statements = self.find_run_statements()
        self.config_files = self.find_config_lines()
        self.nextflow_path = self.find_script_line()
        self.nextflow = NextflowFile(self.nextflow_path)
        self.root_path = self.find_common_path()
        self.dependencies = self.nextflow.includes + self.run_statements
        self.tags = self.find_tags()

    def find_script_line(self) -> Path:
        """
        Find the first script line in the nf-test file.

        Returns:
            Path: The path to the script line.
        """
        for line in self.lines:
            if line.strip().startswith("script"):
                script_path = Path(line.strip().split()[1].strip("\"'"))
                nf_path = self.test_path.parent.joinpath(script_path)
                # If using a relative path
                if nf_path.exists():
                    return nf_path
                # If using relative to the root path
                elif self.repo.joinpath(script_path).exists():
                    return self.repo.joinpath(script_path)
                # Finally, relative to where we're running the script.
                # Unlikely but not impossible
                elif script_path.exists():
                    return script_path
                else:
                    raise FileNotFoundError(
                        f"Script file not found at {nf_path} or {script_path}"
                    )
        else:
            raise ValueError("Script line not found in nf-test file.")

    def find_config_lines(self) -> list[Path]:
        """
        Finds the configuration files mentioned in the lines of the object.

        Returns:
            A list of Path objects representing the paths of the configuration files.
        """
        config_files = []
        for line in self.lines:
            if line.strip().startswith("config"):
                script_path = Path(line.strip().split()[1].strip("\"'"))
                config_path = self.test_path.parent.joinpath(script_path)
                if config_path.exists():
                    config_files.append(config_path.resolve())

        return config_files

    def find_test_type(self) -> tuple[TestTargetType, str]:
        """
        Finds the test type keyword and name from a list of lines.

        Args:
            lines (list[str]): The list of lines to search for the test type.

        Returns:
            tuple(str, str): A tuple containing the test type keyword and name.

        """
        for line in self.lines:
            words = line.split()
            if (
                line.strip().startswith(("workflow", "process", "function"))
                and len(words) == 2
            ):
                keyword = words[0]
                name = words[1].strip("'\"")  # Strip both single and double quotes
                return (TestTargetType(keyword), name)
        return (TestTargetType("pipeline"), "PIPELINE")

    def find_run_statements(self) -> list[str]:
        """
        Find all run statements in a list of lines.

        Args:
            lines (list): List of lines to scan.

        Returns:
            list: List of run statements.
        """
        result = []
        for line in self.lines:
            if line.strip().startswith("run"):
                # This parses `run("<tool>")` to `<tool>
                dependency = (
                    line.strip()
                    .split()[0]
                    .lstrip("run(")
                    .rstrip(")")
                    .strip("\"'")
                    .casefold()
                )
                result.append(dependency)

        return result

    def find_common_path(self) -> Path:
        """
        Finds the common path between the nextflow path and the test path.

        Returns:
            pathlib.Path: The common path between the two paths.
        """
        # Get normalised containing directories
        nf_path_dir = self.nextflow.path.parent.resolve()
        test_path_dir = self.test_path.parent.resolve()
        # Find diff between them
        try:
            diff = nf_path_dir.relative_to(test_path_dir)
        except ValueError:
            diff = nf_path_dir
        # Return common path
        return test_path_dir.joinpath(diff).resolve()

    def find_tags(self) -> list[str]:
        """
        Find all tags in the Nextflow file.
        """
        result = []
        for line in self.lines:
            if line.strip().startswith("tag "):
                logging.debug(f"Found tag line: {line}")
                tag = line.strip().removeprefix("tag ").strip().strip("'\"").casefold()
                result.append(tag)
        logging.debug(f"Found tags: {result}")
        return result

    def detect_if_path_is_in_test(self, path: Path) -> bool:
        """
        Detects if a path is in the test, i.e. the path is either the test itself, the nextflow script, the test directory, or the config file.

        Args:
            path (Path): The path to test if nf-test files match.

        Returns:
            bool: True if the path is in the test files, False otherwise.
        """
        glob_path = self.root_path.resolve().joinpath("*")

        in_root_dir = path.match(glob_path) or path.match(self.root_path.resolve())
        if in_root_dir:
            logging.debug(f"Found {path} in root directory of {self.root_path}")
        match_nf_file = path.match(self.nextflow.path.resolve())
        if match_nf_file:
            logging.debug(f"Found {path} in nextflow file {self.nextflow.path}")
        match_test_file = path.match(self.test_path.resolve())
        if match_test_file:
            logging.debug(f"Found {path} in test file {self.test_path}")
        match_config_file = any(
            path.match(config_file.resolve()) for config_file in self.config_files
        )
        if match_config_file:
            logging.debug(f"Found ${path} in config file {self.config_files}")

        # On the off chance the path is a directory, we need check the files are not within the changed_file path
        in_changed_file_path = any(
            path in x.parents or path == x
            for x in [self.root_path, self.nextflow.path, self.test_path]
            + self.config_files
        )
        if in_changed_file_path:
            logging.debug(f"Found ${path} in test files of {self.test_name}")
        return any(
            [
                in_root_dir,
                match_nf_file,
                match_test_file,
                match_config_file,
                in_changed_file_path,
            ]
        )

    def find_matching_dependencies(self, other_nf_tests):
        """
        Finds and returns a list of NF tests from `other_nf_tests` that have matching test names in `self.dependencies`.

        Args:
            other_nf_tests (list): A list of NF tests to compare against.

        Returns:
            list: A list of NF tests that have matching test names in `self.dependencies`.
        """
        return [
            nf_test
            for nf_test in other_nf_tests
            if nf_test.test_name.casefold() in self.dependencies
        ]

    def get_parents(self, n: int) -> Path:
        """
        Get the parent directory of a path n levels up.
        Args:
            n (int): The number of levels to go up.
        Returns:
            Path: The parent directory n levels up.
        """
        _path = self.test_path
        for _ in range(n):
            _path = _path.parent
        return _path


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
        "-p",
        "--path",
        help="Path to scan for nf-test files. Should be root of repository.",
        default=".",
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
        type=str,
        default="INFO",
        help="Logging level",
    )
    parser.add_argument(
        "-t",
        "--types",
        type=str,
        default="function,process,workflow,pipeline",
        help="Types of tests to include.",
    )
    parser.add_argument(
        "-n",
        "--n_parents",
        type=int,
        default=0,
        help="Number of parents to up to return. 0 for file, 1 for immediate dir, 2 for parent dir, etc.",
    )
    parser.add_argument(
        "-T",
        "--tags",
        type=str,
        default="",
        help="Tags to include.",
    )
    parser.add_argument("--exclude_tags", type=str, default="", help="Tags to exclude.")
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
    path: Path,
    branch1: str,
    branch2: str,
    ignore: list[str],
) -> list[Path]:
    """
    Find all *.nf.tests that are associated with files that have been changed between two specified branches.

    Args:
        repo (Path)        : Path to the repository to scan.
        branch1 (str)      : The first branch being compared
        branch2 (str)      : The second branch being compared
        ignore  (list)     : List of files or file substrings to ignore.

    Returns:
        list: List of files matching the pattern *.nf.test that have changed between branch2 and branch1.
    """

    # Initialise repo for scanning
    repo = Repo(path)

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
        logging.debug(f"File found in diff_index: {str(filepath)}")

        # If file does not match any in the ignore list, add containing directory to changed_files
        match_against_ignored = []
        for ignored_path in ignore:
            is_matched = fnmatch.fnmatch(str(filepath), ignored_path)
            logging.debug(
                f"Checking match of {str(filepath)} against {str(ignored_path)}: {str(is_matched)}"
            )
            match_against_ignored.append(is_matched)
        if not any(match_against_ignored):
            # Prepend the root of the path for better scanning
            changed_files.append(path.joinpath(filepath).resolve())

    # Uniqueify the results before returning for efficiency
    return list(set(changed_files))


def detect_include_files(
    changed_files: list[Path], include_files: dict[str, str], root: Path
) -> list[Path]:
    """
    Detects the include files based on the changed files.

    Args:
        changed_files (list[Path]): List of paths to the changed files.
        include_files (dict[str, str]): Key-value pairs to return if a certain file has changed. If a file in a directory has changed, it points to a different directory.
        root (Path): The root path of the repository.

    Returns:
        list[Path]: List of paths to representing the keys of the include_files dictionary, where a value matched a path in changed_files.
    """
    new_changed_files = []
    for filepath in changed_files:
        # If file is in the include_files, we return the key instead of the value
        for include_path, include_key in include_files.items():
            if filepath.match(include_path):
                new_changed_files.append(root.joinpath(Path(include_key)).resolve())
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


if __name__ == "__main__":
    # Utility stuff
    args = parse_args()
    logging.basicConfig(level=args.log_level)
    # Argparse handling of nargs is a bit rubbish. So we do it manually here.
    args.types = args.types.split(",")

    if args.tags.strip() != "":
        args.tags = [tag.strip().casefold() for tag in args.tags.split(",")]

    if args.exclude_tags.strip() != "":
        args.exclude_tags = [
            tag.strip().casefold() for tag in args.exclude_tags.split(",")
        ]

    # Quick validation of args.types since we cant do this in argparse
    if any(
        _type not in ["function", "process", "workflow", "pipeline"]
        for _type in args.types
    ):
        raise ValueError(
            f"Invalid test type specified. Must be one of 'function', 'process', 'workflow', 'pipeline'. Found: {args.types}"
        )

    root_path = Path(args.path)
    modules_path = Path(root_path, "modules")
    subworkflows_path = Path(root_path, "subworkflows")
    workflows_path = Path(root_path, "workflows")

    logging.info(
        f"Getting files that are different between {args.head_ref} and {args.base_ref}"
    )
    changed_files = find_changed_files(
        root_path, args.head_ref, args.base_ref, args.ignored_files
    )
    logging.debug(f"Found changed files:{changed_files}")

    # If an additional include YAML is added, we detect additional changed dirs to include
    if args.include:
        logging.debug(f"Reading include file: {args.include}")
        include_files = read_yaml_inverted(args.include)
        logging.info(f"Supplementing changed files with include files: {include_files}")
        changed_files = changed_files + detect_include_files(
            changed_files, include_files, root_path
        )

    logging.info("Parsing nf-test files...")
    nf_test_objects = [
        NfTest(_nf_test_file, repo=root_path)
        for _nf_test_file in detect_files([root_path], "*.nf.test")
    ]
    logging.debug(f"Found nf-test files: {[str(x.test_path) for x in nf_test_objects]}")

    logging.info(
        "Getting intersect of changed files and Nextflow components with nf-test files..."
    )
    # Get intersect of changed files and Nextflow components with nf-test files
    logging.info("Finding Nextflow components which have been modified...")
    # changed_files = [changed_files[1]]  # sneaky debugging thing do not merge
    directly_modified_nf_tests = [
        nf_test_object
        for changed_file in changed_files
        for nf_test_object in nf_test_objects
        if nf_test_object.detect_if_path_is_in_test(changed_file)
    ]
    logging.debug(
        f"nf-tests with directly modified nf-test scripts, nextflow scripts, or config files: {[str(x.test_path) for x in directly_modified_nf_tests]}"
    )

    # Get all tests whose dependencies have changed
    logging.info("Finding Nextflow components whose dependencies have changed...")
    indirectly_modified_nf_tests = [
        _nf_test_obj
        for _nf_test_obj in nf_test_objects
        if _nf_test_obj.find_matching_dependencies(directly_modified_nf_tests) != []
    ]
    logging.debug(
        f"Indirectly modified nf-tests: {[str(x.test_path) for x in indirectly_modified_nf_tests]}"
    )

    # Get union of all test files
    logging.debug("Getting union of all nf-test files...")
    all_nf_tests = list(
        {
            nf_test
            for nf_test in directly_modified_nf_tests + indirectly_modified_nf_tests
        }
    )

    # Filter down to only relevant tests
    logging.debug(f"Filtering down to only relevant test types: {args.types}")
    only_selected_nf_tests = [
        nf_test for nf_test in all_nf_tests if nf_test.test_type.value in args.types
    ]
    if args.tags:
        logging.debug(f"Filtering down to only included test tags: {args.tags}")
        logging.debug(
            f"Tests before the filter are: {[test.test_path for test in only_selected_nf_tests]}"
        )
        only_selected_nf_tests = [
            nf_test
            for nf_test in only_selected_nf_tests
            if any(tag in args.tags for tag in nf_test.tags)
        ]
        logging.debug(
            f"Tests after the filter are: {[test.test_path for test in only_selected_nf_tests]}"
        )

    if args.exclude_tags:
        logging.debug(f"Excluding test tags: {args.exclude_tags}")
        logging.debug(
            f"Tests before the filter are: {[test.test_path for test in only_selected_nf_tests]}"
        )
        only_selected_nf_tests = [
            nf_test
            for nf_test in only_selected_nf_tests
            if all(tag not in args.exclude_tags for tag in nf_test.tags)
        ]
        logging.debug(
            f"Tests after the filter are: {[test.test_path for test in only_selected_nf_tests]}"
        )
    # Go back n_parents directories, remove root from path and stringify
    # It's a bit much but might as well do all path manipulation in one place
    logging.info("Normalising test file paths")
    normalised_nf_test_path = sorted(
        list(
            {
                str(
                    nf_test.get_parents(args.n_parents)
                    .resolve()
                    .relative_to(root_path.resolve())
                )
                for nf_test in only_selected_nf_tests
            }
        )
    )

    # Print to string for outputs
    logging.debug("Creating output string...")
    output_string = json.dumps(normalised_nf_test_path)

    if "GITHUB_OUTPUT" in os.environ:
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            print(
                f"components={output_string}",
                file=f,
            )
    logging.info(f"Complete, reporting the following files: {output_string}")
    print(output_string)
