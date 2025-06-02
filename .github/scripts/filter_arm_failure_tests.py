#!/usr/bin/env uv run
# /// script
# dependencies = [
#   "pyyaml",
# ]
# ///

"""
Script to filter nf-test paths based on ARM build failure markers.
This script checks for .arm-build-failure.yml files and filters out tests accordingly.
"""

import json
import yaml
import os
import sys
from pathlib import Path


def has_arm_failure_marker(test_path):
    """
    Check if a test path has an ARM failure marker in its module directory.

    Args:
        test_path: Path to test file or module directory

    Returns:
        tuple: (has_marker, module_dir, failure_info)
    """
    # Extract module directory from test path
    if test_path.startswith("modules/nf-core/"):
        # For module tests like "modules/nf-core/fastqc"
        module_dir = test_path
    elif test_path.startswith("subworkflows/nf-core/"):
        # For subworkflow tests like "subworkflows/nf-core/fastq_align_bwa"
        module_dir = test_path
    else:
        # Try to extract module directory from path
        parts = test_path.split("/")
        if len(parts) >= 3 and parts[0] == "modules" and parts[1] == "nf-core":
            module_dir = "/".join(parts[:3])
        elif len(parts) >= 3 and parts[0] == "subworkflows" and parts[1] == "nf-core":
            module_dir = "/".join(parts[:3])
        else:
            # Can't determine module directory
            return False, None, None

    arm_failure_file = os.path.join(module_dir, ".arm-build-failure.yml")

    if os.path.exists(arm_failure_file):
        try:
            with open(arm_failure_file, 'r') as f:
                failure_info = yaml.safe_load(f) or {}
            return True, module_dir, failure_info
        except Exception as e:
            print(f"Warning: Could not read {arm_failure_file}: {e}", file=sys.stderr)
            return True, module_dir, {}

    return False, module_dir, None


def filter_arm_failure_tests(test_paths):
    """
    Filter test paths based on ARM failure markers.

    Args:
        test_paths: List of test paths

    Returns:
        List of filtered test paths (excluding those with ARM failures)
    """
    filtered_paths = []

    for test_path in test_paths:
        has_failure, module_dir, failure_info = has_arm_failure_marker(test_path)

        if has_failure:
            print(f"⚠️  Skipping tests for {test_path} (ARM build failure detected)", file=sys.stderr)
            if failure_info:
                print(f"   Module: {module_dir}", file=sys.stderr)
                print(f"   Reason: {failure_info.get('reason', 'No reason specified')}", file=sys.stderr)
                print(f"   Date: {failure_info.get('date', 'No date specified')}", file=sys.stderr)
        else:
            filtered_paths.append(test_path)

    return filtered_paths


def main():
    """Main function to handle command line arguments and filter tests."""
    if len(sys.argv) != 2:
        print("Usage: filter_arm_failure_tests.py <test_paths_json>")
        print("test_paths_json: JSON string of test paths")
        sys.exit(1)

    test_paths_json = sys.argv[1]

    try:
        test_paths = json.loads(test_paths_json)
    except json.JSONDecodeError as e:
        print(f"Error parsing test paths JSON: {e}")
        sys.exit(1)

    if not test_paths:
        # Return empty list if no paths
        print(json.dumps([]))
        return

    # Filter tests and output result
    filtered_paths = filter_arm_failure_tests(test_paths)

    print(f"Original test count: {len(test_paths)}", file=sys.stderr)
    print(f"Filtered test count: {len(filtered_paths)}", file=sys.stderr)

    if len(test_paths) != len(filtered_paths):
        skipped_count = len(test_paths) - len(filtered_paths)
        print(f"Skipped {skipped_count} test(s) due to ARM build failures", file=sys.stderr)

    # Output the filtered paths as JSON
    print(json.dumps(filtered_paths))


if __name__ == "__main__":
    main()
