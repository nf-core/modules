import shutil
from unittest import TestCase

import pytest

from nf_core.utils import Pipeline

from .utils import create_tmp_pipeline


class TestPipelines(TestCase):
    def setUp(self) -> None:
        """Create a new Pipeline for testing"""
        self.tmp_dir, self.template_dir, self.pipeline_name, self.pipeline_dir = create_tmp_pipeline()
        self.pipeline_obj = Pipeline(self.pipeline_dir)
        self.pipeline_obj._load()

    def tearDown(self) -> None:
        """Remove the test pipeline directory"""
        shutil.rmtree(self.tmp_dir)

    def _make_pipeline_copy(self):
        """Make a copy of the test pipeline that can be edited

        Returns: Path to new temp directory with pipeline"""
        new_pipeline = self.tmp_dir / "nf-core-testpipeline-copy"
        shutil.copytree(self.pipeline_dir, new_pipeline)
        return new_pipeline

    @pytest.fixture(autouse=True)
    def _use_caplog(self, caplog):
        self.caplog = caplog
