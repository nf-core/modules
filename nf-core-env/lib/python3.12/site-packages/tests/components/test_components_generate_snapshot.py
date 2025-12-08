"""Test generate a snapshot"""

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from nf_core.components.components_test import ComponentsTest
from nf_core.utils import set_wd

from ..test_components import TestComponents
from ..utils import GITLAB_NFTEST_BRANCH, GITLAB_URL


class TestTestComponentsUtils(TestComponents):
    def test_unstable_snapshot(self):
        """Generate the snapshot for a module in nf-core/modules clone with unstable snapshots"""
        with set_wd(self.nfcore_modules):
            snap_generator = ComponentsTest(
                component_type="modules",
                component_name="kallisto/quant",
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
            )
            with pytest.raises(UserWarning) as e:
                snap_generator.run()
            assert "nf-test snapshot is not stable" in str(e.value)

    def test_generate_snapshot_module(self):
        """Generate the snapshot for a module in nf-core/modules clone"""
        with set_wd(self.nfcore_modules):
            snap_generator = ComponentsTest(
                component_type="modules",
                component_name="fastqc",
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
            )
            snap_generator.run()

            snap_path = Path("modules", "nf-core-test", "fastqc", "tests", "main.nf.test.snap")
            assert snap_path.exists()

            with open(snap_path) as fh:
                snap_content = json.load(fh)
            assert "versions" in snap_content
            assert "content" in snap_content["versions"]
            assert "versions.yml:md5,e1cc25ca8af856014824abd842e93978" in snap_content["versions"]["content"][0]

    def test_generate_snapshot_subworkflow(self):
        """Generate the snapshot for a subworkflows in nf-core/modules clone"""
        with set_wd(self.nfcore_modules):
            snap_generator = ComponentsTest(
                component_type="subworkflows",
                component_name="bam_sort_stats_samtools",
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
            )
            snap_generator.run()

            snap_path = Path("subworkflows", "nf-core-test", "bam_sort_stats_samtools", "tests", "main.nf.test.snap")
            assert snap_path.exists()

            with open(snap_path) as fh:
                snap_content = json.load(fh)
            assert "test_bam_sort_stats_samtools_paired_end_flagstats" in snap_content
            assert (
                "test.flagstat:md5,4f7ffd1e6a5e85524d443209ac97d783"
                in snap_content["test_bam_sort_stats_samtools_paired_end_flagstats"]["content"][0][0]
            )
            assert "test_bam_sort_stats_samtools_paired_end_idxstats" in snap_content
            assert (
                "test.idxstats:md5,df60a8c8d6621100d05178c93fb053a2"
                in snap_content["test_bam_sort_stats_samtools_paired_end_idxstats"]["content"][0][0]
            )

    def test_generate_snapshot_once(
        self,
    ):
        """Generate the snapshot for a module in nf-core/modules clone only once"""
        with set_wd(self.nfcore_modules):
            snap_generator = ComponentsTest(
                component_type="modules",
                component_name="fastqc",
                once=True,
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
            )
            snap_generator.repo_type = "modules"
            snap_generator.generate_snapshot = MagicMock()
            snap_generator.run()
            snap_generator.generate_snapshot.assert_called_once()

    def test_update_snapshot_module(self):
        """Update the snapshot of a module in nf-core/modules clone"""

        with set_wd(self.nfcore_modules):
            snap_path = Path("modules", "nf-core-test", "bwa", "mem", "tests", "main.nf.test.snap")
            with open(snap_path) as fh:
                snap_content = json.load(fh)
            original_timestamp = snap_content["Single-End"]["timestamp"]
            # delete the timestamp in json
            snap_content["Single-End"]["timestamp"] = ""
            with open(snap_path, "w") as fh:
                json.dump(snap_content, fh)
            snap_generator = ComponentsTest(
                component_type="modules",
                component_name="bwa/mem",
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
                update=True,
            )
            snap_generator.run()

            with open(snap_path) as fh:
                snap_content = json.load(fh)
            assert "Single-End" in snap_content
            assert snap_content["Single-End"]["timestamp"] != original_timestamp

    def test_test_not_found(self):
        """Generate the snapshot for a module in nf-core/modules clone which doesn't contain tests"""
        with set_wd(self.nfcore_modules):
            snap_generator = ComponentsTest(
                component_type="modules",
                component_name="fastp",
                no_prompts=True,
                remote_url=GITLAB_URL,
                branch=GITLAB_NFTEST_BRANCH,
            )
            test_file = Path("modules", "nf-core-test", "fastp", "tests", "main.nf.test")
            test_file.rename(test_file.parent / "main.nf.test.bak")
            with pytest.raises(UserWarning) as e:
                snap_generator.run()
            assert "Nothing to execute. Is the file 'main.nf.test' missing?" in str(e.value)
            Path(test_file.parent / "main.nf.test.bak").rename(test_file)
