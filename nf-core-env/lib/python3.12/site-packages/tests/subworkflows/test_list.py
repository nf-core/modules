from rich.console import Console

import nf_core.subworkflows

from ..test_subworkflows import TestSubworkflows
from ..utils import GITLAB_SUBWORKFLOWS_BRANCH, GITLAB_URL


class TestSubworkflowsList(TestSubworkflows):
    def test_subworkflows_list_remote(self):
        """Test listing available subworkflows"""
        subworkflows_list = nf_core.subworkflows.SubworkflowList(remote=True)
        listed_subworkflows = subworkflows_list.list_components()
        console = Console(record=True)
        console.print(listed_subworkflows)
        output = console.export_text()
        assert "bam_stats" in output

    def test_subworkflows_list_remote_gitlab(self):
        """Test listing the subworkflows in the remote gitlab repo"""
        subworkflows_list = nf_core.subworkflows.SubworkflowList(
            remote=True, remote_url=GITLAB_URL, branch=GITLAB_SUBWORKFLOWS_BRANCH
        )
        listed_subworkflows = subworkflows_list.list_components()
        console = Console(record=True)
        console.print(listed_subworkflows)
        output = console.export_text()
        assert "bam_stats" in output

    def test_subworkflows_install_and_list_subworkflows(self):
        """Test listing locally installed subworkflows"""
        self.subworkflow_install.install("bam_sort_stats_samtools")
        subworkflows_list = nf_core.subworkflows.SubworkflowList(self.pipeline_dir, remote=False)
        listed_subworkflows = subworkflows_list.list_components()
        console = Console(record=True)
        console.print(listed_subworkflows)
        output = console.export_text()
        assert "bam_stats" in output

    def test_subworkflows_install_gitlab_and_list_subworkflows(self):
        """Test listing locally installed subworkflows"""
        self.subworkflow_install_gitlab.install("bam_sort_stats_samtools")
        subworkflows_list = nf_core.subworkflows.SubworkflowList(self.pipeline_dir, remote=False)
        listed_subworkflows = subworkflows_list.list_components()
        console = Console(record=True)
        console.print(listed_subworkflows)
        output = console.export_text()
        assert "bam_stats" in output
