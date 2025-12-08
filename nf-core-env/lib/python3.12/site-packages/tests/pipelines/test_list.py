"""Tests covering the workflow listing code."""

import json
import os
import tempfile
import time
from datetime import datetime
from pathlib import Path
from unittest import TestCase, mock

import pytest
from rich.console import Console

import nf_core.pipelines.list


class TestList(TestCase):
    """Class for list tests"""

    def setUp(self) -> None:
        # create a temporary directory that can be used by the tests in this file
        tmp = Path(tempfile.TemporaryDirectory().name)
        self.tmp_nxf = tmp / "nxf"
        self.tmp_nxf_str = str(self.tmp_nxf)
        os.environ["NXF_ASSETS"] = self.tmp_nxf_str

    @mock.patch("subprocess.check_output")
    def test_working_listcall(self, mock_subprocess):
        """Test that listing pipelines works"""
        wf_table = nf_core.pipelines.list.list_workflows()
        console = Console(record=True)
        console.print(wf_table)
        output = console.export_text()
        assert "rnaseq" in output
        assert "exoseq" not in output

    @mock.patch("subprocess.check_output")
    def test_working_listcall_archived(self, mock_subprocess):
        """Test that listing pipelines works, showing archived pipelines"""
        wf_table = nf_core.pipelines.list.list_workflows(show_archived=True)
        console = Console(record=True)
        console.print(wf_table)
        output = console.export_text()
        assert "exoseq" in output

    @mock.patch("subprocess.check_output")
    def test_working_listcall_json(self, mock_subprocess):
        """Test that listing pipelines with JSON works"""
        wf_json_str = nf_core.pipelines.list.list_workflows(as_json=True)
        wf_json = json.loads(wf_json_str)
        for wf in wf_json["remote_workflows"]:
            if wf["name"] == "ampliseq":
                break
        else:
            raise AssertionError("Could not find ampliseq in JSON")

    def test_pretty_datetime(self):
        """Test that the pretty datetime function works"""
        now = datetime.now()
        nf_core.pipelines.list.pretty_date(now)
        now_ts = time.mktime(now.timetuple())
        nf_core.pipelines.list.pretty_date(now_ts)

    def test_local_workflows_and_fail(self):
        """Test the local workflow class and try to get local
        Nextflow workflow information"""
        loc_wf = nf_core.pipelines.list.LocalWorkflow("myWF")
        with pytest.raises(RuntimeError):
            loc_wf.get_local_nf_workflow_details()

    def test_local_workflows_compare_and_fail_silently(self):
        """Test the workflow class and try to compare local
        and remote workflows"""
        wfs = nf_core.pipelines.list.Workflows()
        lwf_ex = nf_core.pipelines.list.LocalWorkflow("myWF")
        lwf_ex.full_name = "my Workflow"
        lwf_ex.commit_sha = "aw3s0meh1sh"

        remote = {
            "name": "myWF",
            "full_name": "my Workflow",
            "description": "...",
            "archived": False,
            "stargazers_count": 42,
            "watchers_count": 6,
            "forks_count": 7,
            "releases": [],
        }

        rwf_ex = nf_core.pipelines.list.RemoteWorkflow(remote)
        rwf_ex.commit_sha = "aw3s0meh1sh"
        rwf_ex.releases = [{"tag_sha": "aw3s0meh1sh"}]

        wfs.local_workflows.append(lwf_ex)
        wfs.remote_workflows.append(rwf_ex)
        wfs.compare_remote_local()

        self.assertEqual(rwf_ex.local_wf, lwf_ex)

        rwf_ex.releases = []
        rwf_ex.releases.append({"tag_sha": "noaw3s0meh1sh"})
        wfs.compare_remote_local()

        rwf_ex.full_name = "your Workflow"
        wfs.compare_remote_local()

        rwf_ex.releases = None

    @mock.patch("nf_core.pipelines.list.LocalWorkflow")
    def test_parse_local_workflow_and_succeed(self, mock_local_wf):
        test_path = self.tmp_nxf / "nf-core"
        if not os.path.isdir(test_path):
            os.makedirs(test_path)
        assert os.environ["NXF_ASSETS"] == self.tmp_nxf_str
        with open(self.tmp_nxf / "nf-core/dummy-wf", "w") as f:
            f.write("dummy")
        workflows_obj = nf_core.pipelines.list.Workflows()
        workflows_obj.get_local_nf_workflows()
        assert len(workflows_obj.local_workflows) == 1

    @mock.patch("nf_core.pipelines.list.LocalWorkflow")
    @mock.patch("subprocess.check_output")
    def test_parse_local_workflow_home(self, mock_local_wf, mock_subprocess):
        test_path = self.tmp_nxf / "nf-core"
        if not os.path.isdir(test_path):
            os.makedirs(test_path)
        assert os.environ["NXF_ASSETS"] == self.tmp_nxf_str
        with open(self.tmp_nxf / "nf-core/dummy-wf", "w") as f:
            f.write("dummy")
        workflows_obj = nf_core.pipelines.list.Workflows()
        workflows_obj.get_local_nf_workflows()

    @mock.patch("os.stat")
    @mock.patch("git.Repo")
    def test_local_workflow_investigation(self, mock_repo, mock_stat):
        local_wf = nf_core.pipelines.list.LocalWorkflow("dummy")
        local_wf.local_path = self.tmp_nxf.parent
        mock_repo.head.commit.hexsha = "h00r4y"
        mock_stat.st_mode = 1
        local_wf.get_local_nf_workflow_details()

    def test_worflow_filter(self):
        workflows_obj = nf_core.pipelines.list.Workflows(["rna", "myWF"])

        remote = {
            "name": "myWF",
            "full_name": "my Workflow",
            "description": "rna",
            "archived": False,
            "stargazers_count": 42,
            "watchers_count": 6,
            "forks_count": 7,
            "releases": [],
        }

        rwf_ex = nf_core.pipelines.list.RemoteWorkflow(remote)
        rwf_ex.commit_sha = "aw3s0meh1sh"
        rwf_ex.releases = [{"tag_sha": "aw3s0meh1sh"}]

        remote2 = {
            "name": "myWF",
            "full_name": "my Workflow",
            "description": "dna",
            "archived": False,
            "stargazers_count": 42,
            "watchers_count": 6,
            "forks_count": 7,
            "releases": [],
        }

        rwf_ex2 = nf_core.pipelines.list.RemoteWorkflow(remote2)
        rwf_ex2.commit_sha = "aw3s0meh1sh"
        rwf_ex2.releases = [{"tag_sha": "aw3s0meh1sh"}]

        workflows_obj.remote_workflows.append(rwf_ex)
        workflows_obj.remote_workflows.append(rwf_ex2)

        assert len(workflows_obj.filtered_workflows()) == 1

    def test_filter_archived_workflows(self):
        """
        Test that archived workflows are not shown by default
        """
        workflows_obj = nf_core.pipelines.list.Workflows()
        remote1 = {"name": "myWF", "full_name": "my Workflow", "archived": True, "releases": []}
        rwf_ex1 = nf_core.pipelines.list.RemoteWorkflow(remote1)
        remote2 = {"name": "myWF", "full_name": "my Workflow", "archived": False, "releases": []}
        rwf_ex2 = nf_core.pipelines.list.RemoteWorkflow(remote2)

        workflows_obj.remote_workflows.append(rwf_ex1)
        workflows_obj.remote_workflows.append(rwf_ex2)

        filtered_workflows = workflows_obj.filtered_workflows()
        expected_workflows = [rwf_ex2]

        assert filtered_workflows == expected_workflows

    def test_show_archived_workflows(self):
        """
        Test that archived workflows can be shown optionally
        """
        workflows_obj = nf_core.pipelines.list.Workflows(show_archived=True)
        remote1 = {"name": "myWF", "full_name": "my Workflow", "archived": True, "releases": []}
        rwf_ex1 = nf_core.pipelines.list.RemoteWorkflow(remote1)
        remote2 = {"name": "myWF", "full_name": "my Workflow", "archived": False, "releases": []}
        rwf_ex2 = nf_core.pipelines.list.RemoteWorkflow(remote2)

        workflows_obj.remote_workflows.append(rwf_ex1)
        workflows_obj.remote_workflows.append(rwf_ex2)

        filtered_workflows = workflows_obj.filtered_workflows()
        expected_workflows = [rwf_ex1, rwf_ex2]

        assert filtered_workflows == expected_workflows
