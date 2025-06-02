#!/usr/bin/env uv run
# /// script
# dependencies = [
#   "pyyaml",
# ]
# ///

"""
Script to filter ARM builds based on failure markers.
This script checks for .arm-build-failure.yml files and filters out ARM builds accordingly.
"""

import json
import yaml
import os
import sys
from pathlib import Path


def filter_arm_builds(files_list, build_type="conda"):
    """
    Filter ARM builds based on failure markers.

    Args:
        files_list: List of file paths (environment.yml or Dockerfile)
        build_type: Either "conda" or "dockerfile"

    Returns:
        Dictionary with filtered matrix configuration
    """
    filtered_combinations = []

    for file_path in files_list:
        # Extract module directory from file path
        module_dir = os.path.dirname(file_path)
        arm_failure_file = os.path.join(module_dir, ".arm-build-failure.yml")

        # Check if ARM failure marker exists
        arm_build_disabled = os.path.exists(arm_failure_file)

        if arm_build_disabled:
            print(f"⚠️  ARM builds disabled for {module_dir} (failure marker found)", file=sys.stderr)
            # Read failure info for logging
            try:
                with open(arm_failure_file, 'r') as f:
                    failure_info = yaml.safe_load(f)
                    if failure_info:
                        print(f"   Reason: {failure_info.get('reason', 'No reason specified')}", file=sys.stderr)
                        print(f"   Date: {failure_info.get('date', 'No date specified')}", file=sys.stderr)
            except Exception as e:
                print(f"   Could not read failure info: {e}", file=sys.stderr)

        # Add combinations based on build type
        if build_type == "conda":
            # Add all combinations for conda builds
            for profile in ["docker", "singularity"]:
                for platform in ["linux/amd64", "linux/arm64"]:
                    # Skip ARM builds if failure marker exists
                    if arm_build_disabled and platform == "linux/arm64":
                        print(f"   Skipping {platform} build for {profile}", file=sys.stderr)
                        continue

                    filtered_combinations.append({
                        "files": file_path,
                        "profile": profile,
                        "platform": platform
                    })

        elif build_type == "dockerfile":
            # Add combinations for dockerfile builds
            for platform in ["linux/amd64", "linux/arm64"]:
                # Skip ARM builds if failure marker exists
                if arm_build_disabled and platform == "linux/arm64":
                    print(f"   Skipping {platform} Dockerfile build", file=sys.stderr)
                    continue

                filtered_combinations.append({
                    "files": file_path,
                    "platform": platform
                })

    return {"include": filtered_combinations}


def main():
    """Main function to handle command line arguments and filter builds."""
    if len(sys.argv) != 3:
        print("Usage: script.py <build_type> <files_json>")
        print("build_type: 'conda' or 'dockerfile'")
        print("files_json: JSON string of file paths")
        sys.exit(1)

    build_type = sys.argv[1]
    files_json = sys.argv[2]

    if build_type not in ["conda", "dockerfile"]:
        print("Error: build_type must be 'conda' or 'dockerfile'")
        sys.exit(1)

    try:
        files_list = json.loads(files_json)
    except json.JSONDecodeError as e:
        print(f"Error parsing files JSON: {e}")
        sys.exit(1)

    if not files_list:
        # Return empty matrix if no files
        print(json.dumps({"include": []}))
        return

        # Filter builds and output result
    filtered_matrix = filter_arm_builds(files_list, build_type)
    print(f"Filtered {build_type} matrix: {json.dumps(filtered_matrix, indent=2)}", file=sys.stderr)

    # Output just the JSON for GitHub Actions (to stdout)
    print(json.dumps(filtered_matrix))


if __name__ == "__main__":
    main()
