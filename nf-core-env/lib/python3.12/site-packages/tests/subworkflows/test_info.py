from rich.console import Console

import nf_core.subworkflows

from ..test_subworkflows import TestSubworkflows
from ..utils import GITLAB_SUBWORKFLOWS_BRANCH, GITLAB_URL


class TestSubworkflowsInfo(TestSubworkflows):
    def test_subworkflows_info_remote(self):
        """Test getting info about a remote subworkflow"""
        mods_info = nf_core.subworkflows.SubworkflowInfo(self.pipeline_dir, "bam_sort_stats_samtools")
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Subworkflow: bam_sort_stats_samtools" in output
        assert "Inputs" in output
        assert "Outputs" in output

    def test_subworkflows_info_remote_gitlab(self):
        """Test getting info about a subworkflow in the remote gitlab repo"""
        mods_info = nf_core.subworkflows.SubworkflowInfo(
            self.pipeline_dir, "bam_sort_stats_samtools", remote_url=GITLAB_URL, branch=GITLAB_SUBWORKFLOWS_BRANCH
        )
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Subworkflow: bam_sort_stats_samtools" in output
        assert "Inputs" in output
        assert "Outputs" in output
        assert "--git-remote" in output

    def test_subworkflows_info_local(self):
        """Test getting info about a locally installed subworkflow"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        mods_info = nf_core.subworkflows.SubworkflowInfo(self.pipeline_dir, "bam_sort_stats_samtools")
        mods_info.local = True
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Subworkflow: bam_sort_stats_samtools" in output
        assert "Inputs" in output
        assert "Outputs" in output

    def test_subworkflows_info_in_modules_repo(self):
        """Test getting info about a locally subworkflow in the modules repo"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        mods_info = nf_core.subworkflows.SubworkflowInfo(self.nfcore_modules, "bam_sort_stats_samtools")
        mods_info.local = True
        mods_info_output = mods_info.get_component_info()
        console = Console(record=True)
        console.print(mods_info_output)
        output = console.export_text()

        assert "Subworkflow: bam_sort_stats_samtools" in output
        assert "Inputs" in output
        assert "Outputs" in output
