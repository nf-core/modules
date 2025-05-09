#!/usr/bin/env python3

import argparse
import json
import os
import yaml
import re
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Update meta.yml files with container information')
    parser.add_argument('--conda-files', type=str, default='[]', help='JSON array of changed conda environment.yml files')
    parser.add_argument('--dockerfile-files', type=str, default='[]', help='JSON array of changed Dockerfile files')
    parser.add_argument('--conda-containers', type=str, default='{}', help='JSON object with conda container information')
    parser.add_argument('--dockerfile-containers', type=str, default='{}', help='JSON object with dockerfile container information')
    return parser.parse_args()

def find_meta_yml(file_path):
    """Find the corresponding meta.yml file for an environment.yml or Dockerfile"""
    dir_path = os.path.dirname(file_path)
    meta_path = os.path.join(dir_path, "meta.yml")
    if os.path.exists(meta_path):
        return meta_path
    return None

def parse_container_info(container_info, file_path, container_type, platform):
    """Extract container details from Wave CLI output"""
    if not container_info:
        print(f"No container info found for {file_path}")
        return None

    # Parse container info from Wave output
    try:
        data = json.loads(container_info)
        if not data:
            return None

        # Extract container details
        # This will need to be adjusted based on actual Wave CLI JSON output format
        container_details = {
            "name": data.get("container", {}).get("uri", ""),
            "build_id": data.get("container", {}).get("build_id", ""),
            "scan_id": data.get("container", {}).get("scan_id", "")
        }

        # If it's a singularity container, add the https URL
        if container_type == "singularity":
            container_details["https"] = data.get("container", {}).get("https_url", "")
            # Convert Docker URI to ORAS URI for Singularity
            if container_details["name"] and not container_details["name"].startswith("oras://"):
                container_details["name"] = f"oras://{container_details['name']}"

        return container_details
    except Exception as e:
        print(f"Error parsing container info: {e}")
        return None

def update_meta_yml(meta_path, container_info, container_type, platform):
    """Update the containers section in meta.yml"""
    try:
        with open(meta_path, 'r') as f:
            meta_data = yaml.safe_load(f)

        # Ensure containers section exists
        if 'containers' not in meta_data:
            meta_data['containers'] = {}

        # Ensure container type section exists
        if container_type not in meta_data['containers']:
            meta_data['containers'][container_type] = {}

        # Update platform-specific container info
        meta_data['containers'][container_type][platform] = container_info

        # Add renovate comment above the name field
        name_comment = "# renovate: datasource=docker depName={} versioning=semver".format(
            container_info["name"].split(":")[0]
        )

        # Write back to file with renovate comment preserved
        with open(meta_path, 'w') as f:
            yaml.dump(meta_data, f, default_flow_style=False, sort_keys=False)

        # Add renovate comment back by reading/modifying the file directly
        with open(meta_path, 'r') as f:
            content = f.read()

        # Find the container name line and add the renovate comment above it
        name_pattern = re.compile(rf"({platform}:(?:\s*)\n(?:\s*)name: {re.escape(container_info['name'])})")
        if name_pattern.search(content):
            content = name_pattern.sub(rf"\1\n            {name_comment}", content)

        with open(meta_path, 'w') as f:
            f.write(content)

        print(f"Updated {meta_path} with {container_type} {platform} container info")
    except Exception as e:
        print(f"Error updating {meta_path}: {e}")

def main():
    args = parse_args()

    # Parse input arguments
    conda_files = json.loads(args.conda_files)
    dockerfile_files = json.loads(args.dockerfile_files)
    conda_containers = args.conda_containers
    dockerfile_containers = args.dockerfile_containers

    # Process conda environment.yml files
    for file in conda_files:
        meta_path = find_meta_yml(file)
        if not meta_path:
            print(f"No meta.yml found for {file}")
            continue

        # Get module name to match with container info
        module_name = os.path.dirname(file).split('/')[-1]

        # Update for Docker and Singularity containers
        for container_type in ["docker", "singularity"]:
            for platform in ["linux_amd64", "linux_arm64"]:
                container_info = parse_container_info(
                    conda_containers,
                    file,
                    container_type,
                    platform
                )
                if container_info:
                    update_meta_yml(meta_path, container_info, container_type, platform)

    # Process Dockerfile files
    for file in dockerfile_files:
        meta_path = find_meta_yml(file)
        if not meta_path:
            print(f"No meta.yml found for {file}")
            continue

        # Get module name to match with container info
        module_name = os.path.dirname(file).split('/')[-1]

        # Update for Docker containers only
        for platform in ["linux_amd64", "linux_arm64"]:
            container_info = parse_container_info(
                dockerfile_containers,
                file,
                "docker",
                platform
            )
            if container_info:
                update_meta_yml(meta_path, container_info, "docker", platform)

if __name__ == "__main__":
    main()
