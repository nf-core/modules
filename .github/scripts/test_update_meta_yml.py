#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.8"
# dependencies = [
#   "pytest>=7.0.0",
#   "pytest-cov>=4.0.0",
#   "pyyaml>=6.0",
# ]
# ///

"""
Self-contained test script for update_meta_yml.py

Usage:
  chmod +x test_update_meta_yml.py
  ./test_update_meta_yml.py
"""

import pytest
import os
import sys
import json
import yaml
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, mock_open, MagicMock

# Add the current directory to the path so we can import the script modules
sys.path.insert(0, str(Path(__file__).parent))

# Import functions from the script to test
from update_meta_yml import find_meta_yml, parse_container_info, update_meta_yml, parse_args, main

# Fixtures
@pytest.fixture
def sample_container_info():
    return json.dumps({
        "container": {
            "uri": "community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42",
            "build_id": "5cfd0f3cb6760c42_1",
            "scan_id": "6fc310277b74",
            "https_url": "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b280a35770a70ed67008c1d6b6db118409bc3adbb3a98edcd55991189e5116f6/data"
        }
    })

@pytest.fixture
def sample_meta_yml_content():
    return {
        "name": "fastqc",
        "description": "Run FastQC on sequenced reads",
        "containers": {
            "docker": {
                "linux_amd64": {
                    "name": "community.wave.seqera.io/library/fastqc:0.12.0--old_hash",
                    "build_id": "old_build_id",
                    "scan_id": "old_scan_id"
                }
            }
        }
    }

@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing file operations"""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

# Tests for find_meta_yml
def test_find_meta_yml_exists(temp_dir):
    # Create a module directory structure
    module_dir = os.path.join(temp_dir, "modules", "nf-core", "fastqc")
    os.makedirs(module_dir, exist_ok=True)

    # Create environment.yml and meta.yml files
    env_file = os.path.join(module_dir, "environment.yml")
    meta_file = os.path.join(module_dir, "meta.yml")

    with open(env_file, 'w') as f:
        f.write("# Test environment file")

    with open(meta_file, 'w') as f:
        f.write("# Test meta file")

    # Test finding meta.yml from environment.yml
    result = find_meta_yml(env_file)
    assert result == meta_file

def test_find_meta_yml_not_exists(temp_dir):
    # Create a module directory without meta.yml
    module_dir = os.path.join(temp_dir, "modules", "nf-core", "fastqc")
    os.makedirs(module_dir, exist_ok=True)

    env_file = os.path.join(module_dir, "environment.yml")
    with open(env_file, 'w') as f:
        f.write("# Test environment file")

    # Test finding meta.yml that doesn't exist
    result = find_meta_yml(env_file)
    assert result is None

# Tests for parse_container_info
def test_parse_container_info_docker(sample_container_info):
    file_path = "modules/nf-core/fastqc/environment.yml"
    container_type = "docker"
    platform = "linux_amd64"

    result = parse_container_info(sample_container_info, file_path, container_type, platform)

    assert result is not None
    assert result["name"] == "community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42"
    assert result["build_id"] == "5cfd0f3cb6760c42_1"
    assert result["scan_id"] == "6fc310277b74"
    assert "https" not in result

def test_parse_container_info_singularity(sample_container_info):
    file_path = "modules/nf-core/fastqc/environment.yml"
    container_type = "singularity"
    platform = "linux_amd64"

    result = parse_container_info(sample_container_info, file_path, container_type, platform)

    assert result is not None
    assert result["name"] == "oras://community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42"
    assert result["build_id"] == "5cfd0f3cb6760c42_1"
    assert result["scan_id"] == "6fc310277b74"
    assert result["https"] == "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b280a35770a70ed67008c1d6b6db118409bc3adbb3a98edcd55991189e5116f6/data"

def test_parse_container_info_invalid_json():
    file_path = "modules/nf-core/fastqc/environment.yml"
    container_type = "docker"
    platform = "linux_amd64"
    invalid_json = "this is not json"

    result = parse_container_info(invalid_json, file_path, container_type, platform)
    assert result is None

def test_parse_container_info_empty():
    file_path = "modules/nf-core/fastqc/environment.yml"
    container_type = "docker"
    platform = "linux_amd64"

    result = parse_container_info("", file_path, container_type, platform)
    assert result is None

# Tests for update_meta_yml
def test_update_meta_yml(temp_dir, sample_meta_yml_content):
    # Create test meta.yml file
    meta_file = os.path.join(temp_dir, "meta.yml")
    with open(meta_file, 'w') as f:
        yaml.dump(sample_meta_yml_content, f)

    # Container info to update
    container_info = {
        "name": "community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42",
        "build_id": "5cfd0f3cb6760c42_1",
        "scan_id": "6fc310277b74"
    }

    # Update meta.yml
    update_meta_yml(meta_file, container_info, "docker", "linux_amd64")

    # Read updated file
    with open(meta_file, 'r') as f:
        updated_content = yaml.safe_load(f)

    # Verify changes
    assert updated_content["containers"]["docker"]["linux_amd64"]["name"] == container_info["name"]
    assert updated_content["containers"]["docker"]["linux_amd64"]["build_id"] == container_info["build_id"]
    assert updated_content["containers"]["docker"]["linux_amd64"]["scan_id"] == container_info["scan_id"]

def test_update_meta_yml_new_container_type(temp_dir, sample_meta_yml_content):
    # Create test meta.yml file
    meta_file = os.path.join(temp_dir, "meta.yml")
    with open(meta_file, 'w') as f:
        yaml.dump(sample_meta_yml_content, f)

    # Container info to update (for a new container type)
    container_info = {
        "name": "oras://community.wave.seqera.io/library/fastqc:0.12.1--0827550dd72a3745",
        "build_id": "0827550dd72a3745_1",
        "scan_id": "some_scan_id",
        "https": "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b280a35770a70ed67008c1d6b6db118409bc3adbb3a98edcd55991189e5116f6/data"
    }

    # Update meta.yml with singularity container
    update_meta_yml(meta_file, container_info, "singularity", "linux_amd64")

    # Read updated file
    with open(meta_file, 'r') as f:
        updated_content = yaml.safe_load(f)

    # Verify changes - new container type added
    assert "singularity" in updated_content["containers"]
    assert updated_content["containers"]["singularity"]["linux_amd64"]["name"] == container_info["name"]
    assert updated_content["containers"]["singularity"]["linux_amd64"]["build_id"] == container_info["build_id"]
    assert updated_content["containers"]["singularity"]["linux_amd64"]["scan_id"] == container_info["scan_id"]
    assert updated_content["containers"]["singularity"]["linux_amd64"]["https"] == container_info["https"]

# Tests for main function using mocks
@patch('update_meta_yml.parse_args')
@patch('update_meta_yml.find_meta_yml')
@patch('update_meta_yml.parse_container_info')
@patch('update_meta_yml.update_meta_yml')
def test_main_conda_files(mock_update_meta, mock_parse_container, mock_find_meta, mock_parse_args):
    # Setup mocks
    mock_args = MagicMock()
    mock_args.conda_files = json.dumps(["modules/nf-core/fastqc/environment.yml"])
    mock_args.dockerfile_files = json.dumps([])
    mock_args.conda_containers = "{}"
    mock_args.dockerfile_containers = "{}"
    mock_parse_args.return_value = mock_args

    mock_find_meta.return_value = "modules/nf-core/fastqc/meta.yml"
    mock_parse_container.return_value = {
        "name": "community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42",
        "build_id": "5cfd0f3cb6760c42_1",
        "scan_id": "6fc310277b74"
    }

    # Call main
    main()

    # Verify the right functions were called with correct arguments
    assert mock_find_meta.call_count == 1
    assert mock_parse_container.call_count == 4  # 2 container types * 2 platforms
    assert mock_update_meta.call_count == 4

# Run the tests when script is executed directly
if __name__ == "__main__":
    # Use sys.argv[1:] to pass any command-line arguments to pytest
    exit_code = pytest.main(["-v", __file__] + sys.argv[1:])
    sys.exit(exit_code)
