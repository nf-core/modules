"""Tests covering the refgenie integration code"""

import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path


class TestRefgenie(unittest.TestCase):
    """Class for refgenie tests"""

    def setUp(self):
        """
        Prepare a refgenie config file
        """
        self.tmp_dir = Path(tempfile.TemporaryDirectory().name)
        self.NXF_HOME = self.tmp_dir / ".nextflow"
        self.NXF_REFGENIE_PATH = self.NXF_HOME / "nf-core" / "refgenie_genomes.config"
        self.REFGENIE = self.tmp_dir / "genomes_config.yaml"
        self.translation_file = self.tmp_dir / "alias_translations.yaml"
        # Set NXF_HOME environment variable
        # avoids adding includeConfig statement to config file outside the current tmpdir
        try:
            self.NXF_HOME_ORIGINAL = os.environ["NXF_HOME"]
        except Exception:
            self.NXF_HOME_ORIGINAL = None
        os.environ["NXF_HOME"] = str(self.NXF_HOME)

        # create NXF_HOME and nf-core directories
        nf_core_dir = self.NXF_HOME / "nf-core"
        nf_core_dir.mkdir(parents=True, exist_ok=True)

        # Initialize a refgenie config
        os.system(f"refgenie init -c {self.REFGENIE}")

        # Add NXF_REFGENIE_PATH to refgenie config
        with open(self.REFGENIE, "a") as fh:
            fh.write(f"nextflow_config: {self.NXF_REFGENIE_PATH}\n")

        # Add an alias translation to YAML file
        with open(self.translation_file, "a") as fh:
            fh.write("ensembl_gtf: gtf\n")

    def tearDown(self) -> None:
        # Reset NXF_HOME environment variable
        if self.NXF_HOME_ORIGINAL is None:
            del os.environ["NXF_HOME"]
        else:
            os.environ["NXF_HOME"] = self.NXF_HOME_ORIGINAL

    def test_update_refgenie_genomes_config(self):
        """Test that listing pipelines works"""
        # Populate the config with a genome
        cmd = f"refgenie pull t7/fasta -c {self.REFGENIE}"
        out = subprocess.check_output(shlex.split(cmd), stderr=subprocess.STDOUT)

        assert "Updated nf-core genomes config" in str(out)

    def test_asset_alias_translation(self):
        """Test that asset aliases are translated correctly"""
        # Populate the config with a genome
        cmd = f"refgenie pull hg38/ensembl_gtf -c {self.REFGENIE}"
        subprocess.check_output(shlex.split(cmd), stderr=subprocess.STDOUT)
        cmd = f"cat {self.NXF_REFGENIE_PATH}"
        out = subprocess.check_output(shlex.split(cmd), stderr=subprocess.STDOUT)
        assert "      gtf                  = " in str(out)
        assert "      ensembl_gtf          = " not in str(out)
