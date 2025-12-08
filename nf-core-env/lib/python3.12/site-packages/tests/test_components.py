"""Tests covering the modules commands"""

import os
import shutil
import tempfile
import unittest
from pathlib import Path

import pytest
from git.repo import Repo

from .utils import GITLAB_NFTEST_BRANCH, GITLAB_URL


class TestComponents(unittest.TestCase):
    """Class for components tests"""

    def setUp(self):
        """Clone a testing version the nf-core/modules repo"""
        self.tmp_dir = Path(tempfile.mkdtemp())
        self.nfcore_modules = Path(self.tmp_dir, "modules-test")

        Repo.clone_from(GITLAB_URL, self.nfcore_modules, branch=GITLAB_NFTEST_BRANCH)

        # Set $PROFILE environment variable to docker - tests will run with Docker
        if os.environ.get("PROFILE") is None:
            os.environ["PROFILE"] = "docker"

    def tearDown(self):
        """Clean up temporary files and folders"""

        # Clean up temporary files
        if self.tmp_dir.is_dir():
            shutil.rmtree(self.tmp_dir, ignore_errors=True)

    @pytest.fixture(autouse=True)
    def _use_caplog(self, caplog):
        self.caplog = caplog
