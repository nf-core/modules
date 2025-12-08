"""Global pytest configuration for nf-core tests setting up worker-specific cache directories to avoid git lock issues."""

import os
import shutil
import tempfile


def pytest_configure(config):
    """Configure pytest before any tests run - set up worker-specific cache directories."""
    # Get worker ID for pytest-xdist, or 'main' if not using xdist
    worker_id = getattr(config, "workerinput", {}).get("workerid", "main")

    # Create temporary directories for this worker
    cache_base = tempfile.mkdtemp(prefix=f"nfcore_cache_{worker_id}_")
    config_base = tempfile.mkdtemp(prefix=f"nfcore_config_{worker_id}_")

    # Store original values for later restoration
    config._original_xdg_cache = os.environ.get("XDG_CACHE_HOME")
    config._original_xdg_config = os.environ.get("XDG_CONFIG_HOME")
    config._temp_cache_dir = cache_base
    config._temp_config_dir = config_base

    # Set environment variables to use worker-specific directories
    os.environ["XDG_CACHE_HOME"] = cache_base
    os.environ["XDG_CONFIG_HOME"] = config_base


def pytest_unconfigure(config):
    """Clean up after all tests are done."""
    # Restore original environment variables
    if hasattr(config, "_original_xdg_cache"):
        if config._original_xdg_cache is not None:
            os.environ["XDG_CACHE_HOME"] = config._original_xdg_cache
        else:
            os.environ.pop("XDG_CACHE_HOME", None)

    if hasattr(config, "_original_xdg_config"):
        if config._original_xdg_config is not None:
            os.environ["XDG_CONFIG_HOME"] = config._original_xdg_config
        else:
            os.environ.pop("XDG_CONFIG_HOME", None)

    # Clean up temporary directories
    if hasattr(config, "_temp_cache_dir"):
        try:
            shutil.rmtree(config._temp_cache_dir)
        except (OSError, FileNotFoundError):
            pass
    if hasattr(config, "_temp_config_dir"):
        try:
            shutil.rmtree(config._temp_config_dir)
        except (OSError, FileNotFoundError):
            pass
