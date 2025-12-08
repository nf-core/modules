"""Tests covering the sync command"""

import json
import os
from pathlib import Path
from unittest import mock

import git
import pytest
import yaml

import nf_core.pipelines.create.create
import nf_core.pipelines.sync
from nf_core.utils import NFCoreYamlConfig

from ..test_pipelines import TestPipelines
from ..utils import with_temporary_folder


class MockResponse:
    def __init__(self, data: dict | list[dict], status_code: int, url: str):
        self.url: str = url
        self.status_code: int = status_code
        self.from_cache: bool = False
        self.reason: str = "Mocked response"
        self.data: dict | list[dict] = data
        self.content: str = json.dumps(data)
        self.headers: dict[str, str] = {"content-encoding": "test", "connection": "fake"}

    def json(self):
        return self.data


def mocked_requests_get(url, params=None, **kwargs) -> MockResponse:
    """Helper function to emulate GET request responses from the web"""

    url_template = "https://api.github.com/repos/{}/response/"
    if url == Path(url_template.format("no_existing_pr"), "pulls?head=TEMPLATE&base=None"):
        return MockResponse([], 200, url)
    if url == Path(url_template.format("list_prs"), "pulls"):
        response_data = [
            {
                "state": "closed",
                "head": {"ref": "nf-core-template-merge-2"},
                "base": {"ref": "main"},
                "html_url": "pr_url",
            }
        ] + [
            {
                "state": "open",
                "head": {"ref": f"nf-core-template-merge-{branch_no}"},
                "base": {"ref": "main"},
                "html_url": "pr_url",
            }
            for branch_no in range(3, 7)
        ]
        return MockResponse(response_data, 200, url)
    if url == "https://nf-co.re/pipelines.json":
        return MockResponse({"remote_workflows": [{"name": "testpipeline", "topics": ["test", "pipeline"]}]}, 200, url)

    return MockResponse([{"html_url": url}], 404, url)


def mocked_requests_patch(url: str, data: str, **kwargs) -> MockResponse:
    """Helper function to emulate POST requests responses from the web"""

    if url == "url_to_update_pr":
        return MockResponse({"html_url": "great_success"}, 200, url)
    # convert data to dict
    response = json.loads(data)
    response["patch_url"] = url
    return MockResponse(response, 404, url)


def mocked_requests_post(url, **kwargs):
    """Helper function to emulate POST requests responses from the web"""

    if url == "https://api.github.com/repos/no_existing_pr/response/pulls":
        return MockResponse({"html_url": "great_success"}, 201, url)

    return MockResponse({}, 404, url)


class TestModules(TestPipelines):
    """Class for modules tests"""

    def setUp(self):
        super().setUp()
        self.remote_path = Path(self.tmp_dir, "remote_repo")
        self.remote_repo = git.Repo.init(self.remote_path, bare=True)

        if self.remote_repo.active_branch.name != "master":
            self.remote_repo.active_branch.rename("master")

    @with_temporary_folder
    def test_inspect_sync_dir_notgit(self, tmp_dir: str):
        """Try syncing an empty directory"""
        nf_core_yml_path = Path(tmp_dir, ".nf-core.yml")
        nf_core_yml = NFCoreYamlConfig(repository_type="pipeline")

        with open(nf_core_yml_path, "w") as fh:
            yaml.dump(nf_core_yml.model_dump(), fh)

        psync = nf_core.pipelines.sync.PipelineSync(tmp_dir)
        with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
            psync.inspect_sync_dir()
        assert "does not appear to be a git repository" in exc_info.value.args[0]

    def test_inspect_sync_dir_dirty(self):
        """Try syncing a pipeline with uncommitted changes"""
        # Add an empty file, uncommitted
        test_fn = Path(self.pipeline_dir) / "uncommitted"
        test_fn.touch()
        # Try to sync, check we halt with the right error
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        try:
            with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
                psync.inspect_sync_dir()
            assert exc_info.value.args[0].startswith("Uncommitted changes found in pipeline directory!")
        finally:
            os.remove(test_fn)

    def test_inspect_sync_ignored_files(self):
        """
        Try inspecting the repo for syncing with untracked changes that are ignored.
        No assertions, we are checking that no exception is raised in the process.
        """
        test_fn = Path(self.pipeline_dir) / "ignored.txt"
        self._make_ignored_file(test_fn)

        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()

    def test_get_wf_config_no_branch(self):
        """Try getting a workflow config when the branch doesn't exist"""
        # Try to sync, check we halt with the right error
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir, from_branch="foo")
        with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
            psync.inspect_sync_dir()
            psync.get_wf_config()
        assert exc_info.value.args[0] == "Branch `foo` not found!"

    def test_get_wf_config_missing_required_config(self):
        """Try getting a workflow config, then make it miss a required config option"""
        # Try to sync, check we halt with the right error
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.required_config_vars = ["fakethisdoesnotexist"]
        with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
            psync.inspect_sync_dir()
            psync.get_wf_config()
        # Check that we did actually get some config back
        assert psync.wf_config["params.validate_params"] == "true"
        # Check that we raised because of the missing fake config var
        assert exc_info.value.args[0] == "Workflow config variable `fakethisdoesnotexist` not found!"

    def test_checkout_template_branch(self):
        """Try checking out the TEMPLATE branch of the pipeline"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()

    def test_checkout_template_branch_no_template(self):
        """Try checking out the TEMPLATE branch of the pipeline when it does not exist"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()

        psync.repo.delete_head("TEMPLATE")

        with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
            psync.checkout_template_branch()
        assert exc_info.value.args[0] == "Could not check out branch 'origin/TEMPLATE' or 'TEMPLATE'"

    def test_delete_tracked_template_branch_files(self):
        """Confirm that we can delete all tracked files in the TEMPLATE branch"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        psync.delete_tracked_template_branch_files()
        top_level_ignored = self._get_top_level_ignored(psync)
        assert set(os.listdir(self.pipeline_dir)) == set([".git"]).union(top_level_ignored)

    def test_delete_tracked_template_branch_files_unlink_throws_error(self):
        """Test that SyncExceptionError is raised when Path.unlink throws an exception"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()

        # Create a test file that would normally be deleted
        test_file = Path(self.pipeline_dir) / "test_file.txt"
        test_file.touch()

        # Mock Path.unlink to raise an exception
        with mock.patch("pathlib.Path.unlink", side_effect=OSError("Permission denied")) as mock_unlink:
            with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
                psync.delete_tracked_template_branch_files()

            # Verify the exception contains the original error
            assert "Permission denied" in str(exc_info.value)

            # Verify Path.unlink was called
            mock_unlink.assert_called()

    def test_delete_tracked_template_branch_rmdir_throws_error(self):
        """Test that SyncExceptionError is raised when Path.rmdir throws an exception"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()

        # Create an empty directory that would normally be deleted
        empty_dir = Path(self.pipeline_dir) / "empty_test_dir"
        empty_dir.mkdir()

        # Mock Path.rmdir to raise an exception
        with mock.patch("pathlib.Path.rmdir", side_effect=OSError("Permission denied")) as mock_rmdir:
            with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
                psync.delete_tracked_template_branch_files()

            # Verify the exception contains the original error
            assert "Permission denied" in str(exc_info.value)

            # Verify Path.rmdir was called
            mock_rmdir.assert_called()

    def test_delete_staged_template_branch_files_ignored(self):
        """Confirm that files in .gitignore are not deleted by delete_staged_template_branch_files"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)

        ignored_file = Path(self.pipeline_dir) / "ignored.txt"
        self._make_ignored_file(ignored_file)

        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        psync.delete_tracked_template_branch_files()

        # Ignored file should still exist
        assert ignored_file.exists()

        # .git directory should still exist
        assert (Path(self.pipeline_dir) / ".git").exists()

    def test_delete_staged_template_branch_files_ignored_nested_dir(self):
        """Confirm that deletion of ignored files respects directory structure"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        repo = git.Repo(self.pipeline_dir)

        # Create this structure:
        # dir
        # ├── subdirA                   # (should be kept)
        # │   └── subdirB               # (should be kept)
        # │       └── ignored.txt       # add to .gitignore (should be kept)
        # └── subdirC                   # (should be deleted)
        #     └── subdirD               # (should be deleted)
        #         └── not_ignored.txt   # commit this file (should be deleted)
        parent_dir = Path(self.pipeline_dir) / "dir"
        to_be_kept_dir = parent_dir / "subdirA" / "subdirB"
        ignored_file = to_be_kept_dir / "ignored.txt"
        to_be_deleted_dir = parent_dir / "subdirC" / "subdirD"
        non_ignored_file = to_be_deleted_dir / "not_ignored.txt"

        to_be_kept_dir.mkdir(parents=True)
        to_be_deleted_dir.mkdir(parents=True)
        non_ignored_file.touch()

        repo.git.add(non_ignored_file)
        repo.index.commit("Add non-ignored file")

        self._make_ignored_file(ignored_file)

        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        psync.delete_tracked_template_branch_files()

        # Ignored file and its parent directory should still exist
        assert ignored_file.exists()
        assert to_be_kept_dir.exists()  # subdirB
        assert to_be_kept_dir.parent.exists()  # subdirA

        # Non-ignored file and its parent directory should be deleted
        assert not non_ignored_file.exists()
        assert not to_be_deleted_dir.exists()  # subdirD
        assert not to_be_deleted_dir.parent.exists()  # subdirC

        # .git directory should still exist
        assert (Path(self.pipeline_dir) / ".git").exists()

    def test_create_template_pipeline(self):
        """Confirm that we can create a new template pipeline in an empty directory"""
        # First, delete all the files
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        psync.delete_tracked_template_branch_files()
        top_level_ignored = self._get_top_level_ignored(psync)
        assert set(os.listdir(self.pipeline_dir)) == set([".git"]).union(top_level_ignored)
        # Now create the new template
        psync.make_template_pipeline()
        assert "main.nf" in os.listdir(self.pipeline_dir)
        assert "nextflow.config" in os.listdir(self.pipeline_dir)

    def test_commit_template_changes_nochanges(self):
        """Try to commit the TEMPLATE branch, but no changes were made"""
        # Check out the TEMPLATE branch but skip making the new template etc.
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        # Function returns False if no changes were made
        assert psync.commit_template_changes() is False

    def test_commit_template_changes_changes(self):
        """Try to commit the TEMPLATE branch, but no changes were made"""
        # Check out the TEMPLATE branch but skip making the new template etc.
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        # Add an empty file, uncommitted
        test_fn = Path(self.pipeline_dir) / "uncommitted"
        test_fn.touch()
        # Check that we have uncommitted changes
        assert psync.repo.is_dirty(untracked_files=True) is True
        # Function returns True if no changes were made
        assert psync.commit_template_changes() is True
        # Check that we don't have any uncommitted changes
        assert psync.repo.is_dirty(untracked_files=True) is False

    def test_commit_template_preserves_ignored(self):
        """Try to commit the TEMPLATE branch, but no changes were made"""
        # Check out the TEMPLATE branch but skip making the new template etc.
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)

        ignored_file = Path(self.pipeline_dir) / "ignored.txt"

        self._make_ignored_file(ignored_file)

        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        psync.commit_template_changes()

        assert ignored_file.exists()

    def test_push_template_branch_error(self):
        """Try pushing the changes, but without a remote (should fail)"""
        # Check out the TEMPLATE branch but skip making the new template etc.
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.checkout_template_branch()
        # Add an empty file and commit it
        test_fn = Path(self.pipeline_dir) / "uncommitted"
        test_fn.touch()
        psync.commit_template_changes()
        # Try to push changes
        with pytest.raises(nf_core.pipelines.sync.PullRequestExceptionError) as exc_info:
            psync.push_template_branch()
        assert exc_info.value.args[0].startswith("Could not push TEMPLATE branch")

    def test_create_merge_base_branch(self):
        """Try creating a merge base branch"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()

        psync.create_merge_base_branch()

        assert psync.merge_branch in psync.repo.branches

    def test_create_merge_base_branch_thrice(self):
        """Try creating a merge base branch thrice

        This is needed because the first time this function is called, the
        merge branch does not exist yet (it is only created at the end of the
        create_merge_base_branch function) and the if-statement is ignored.
        Also, the second time this function is called, the existing merge
        branch only has the base format, i.e. without the -{branch_no} at the
        end, so it is needed to call it a third time to make sure this is
        picked up.
        """
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()

        for _ in range(3):
            psync.create_merge_base_branch()

        assert psync.merge_branch in psync.repo.branches
        for branch_no in [2, 3]:
            assert f"{psync.original_merge_branch}-{branch_no}" in psync.repo.branches

    def test_push_merge_branch(self):
        """Try pushing merge branch"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.repo.create_remote("origin", self.remote_path)

        psync.create_merge_base_branch()
        psync.push_merge_branch()

        assert psync.merge_branch in [b.name for b in self.remote_repo.branches]

    def test_push_merge_branch_without_create_branch(self):
        """Try pushing merge branch without creating first"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()
        psync.repo.create_remote("origin", self.remote_path)

        with pytest.raises(nf_core.pipelines.sync.PullRequestExceptionError) as exc_info:
            psync.push_merge_branch()
        assert exc_info.value.args[0].startswith(f"Could not push branch '{psync.merge_branch}'")

    @mock.patch("nf_core.utils.gh_api.get", side_effect=mocked_requests_get)
    @mock.patch("nf_core.utils.gh_api.post", side_effect=mocked_requests_post)
    def test_make_pull_request_success(self, mock_post, mock_get):
        """Try making a PR - successful response"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.gh_api.get = mock_get
        psync.gh_api.post = mock_post
        psync.gh_username = "no_existing_pr"
        psync.gh_repo = "no_existing_pr/response"
        os.environ["GITHUB_AUTH_TOKEN"] = "test"
        psync.make_pull_request()
        assert psync.gh_pr_returned_data["html_url"] == "great_success"

    @mock.patch("nf_core.utils.gh_api.get", side_effect=mocked_requests_get)
    @mock.patch("nf_core.utils.gh_api.post", side_effect=mocked_requests_post)
    def test_make_pull_request_bad_response(self, mock_post, mock_get):
        """Try making a PR and getting a 404 error"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.gh_api.get = mock_get
        psync.gh_api.post = mock_post
        psync.gh_username = "bad_url"
        psync.gh_repo = "bad_url/response"
        os.environ["GITHUB_AUTH_TOKEN"] = "test"
        with pytest.raises(nf_core.pipelines.sync.PullRequestExceptionError) as exc_info:
            psync.make_pull_request()
        assert exc_info.value.args[0].startswith(
            "Something went badly wrong - GitHub API PR failed - got return code 404"
        )

    def test_reset_target_dir(self):
        """Try resetting target pipeline directory"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()

        psync.repo.git.checkout("dev")

        psync.reset_target_dir()

        assert psync.repo.heads[0].name == "TEMPLATE"

    def test_reset_target_dir_fake_branch(self):
        """Try resetting target pipeline directory but original branch does not exist"""
        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)
        psync.inspect_sync_dir()
        psync.get_wf_config()

        psync.original_branch = "fake_branch"

        with pytest.raises(nf_core.pipelines.sync.SyncExceptionError) as exc_info:
            psync.reset_target_dir()
        assert exc_info.value.args[0].startswith("Could not reset to original branch `fake_branch`")

    def test_sync_no_changes(self):
        """Test pipeline sync when no changes are needed"""
        with (
            mock.patch("requests.get", side_effect=mocked_requests_get),
            mock.patch("requests.post", side_effect=mocked_requests_post) as mock_post,
        ):
            psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)

            # Mock that no changes were made
            psync.made_changes = False

            # Run sync
            psync.sync()

            # Verify no PR was created
            mock_post.assert_not_called()

    def test_sync_no_github_token(self):
        """Test sync fails appropriately when GitHub token is missing"""
        # Ensure GitHub token is not set
        if "GITHUB_AUTH_TOKEN" in os.environ:
            del os.environ["GITHUB_AUTH_TOKEN"]

        psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir, make_pr=True)
        psync.made_changes = True  # Force changes to trigger PR attempt

        # Run sync and check for appropriate error
        with self.assertRaises(nf_core.pipelines.sync.PullRequestExceptionError) as exc_info:
            psync.sync()
        self.assertIn("GITHUB_AUTH_TOKEN not set!", str(exc_info.exception))

    def test_sync_preserves_ignored_files(self):
        """Test that sync preserves files and directories specified in .gitignore"""
        with (
            mock.patch("requests.get", side_effect=mocked_requests_get),
        ):
            psync = nf_core.pipelines.sync.PipelineSync(self.pipeline_dir)

            ignored_file = Path(self.pipeline_dir) / "ignored.txt"
            self._make_ignored_file(ignored_file)

            psync.made_changes = True

            psync.sync()

            self.assertTrue(ignored_file.exists())

    def _make_ignored_file(self, file_path: Path):
        """Helper function to create an ignored file."""
        if not self.pipeline_dir:
            raise ValueError("Instantiate a pipeline before adding ignored files.")

        file_path.touch()

        gitignore_path = Path(self.pipeline_dir) / ".gitignore"
        with open(gitignore_path, "a") as f:
            f.write(f"{file_path.name}\n")

        repo = git.Repo(self.pipeline_dir)
        repo.git.add(".gitignore")
        repo.index.commit("Add .gitignore")

    def _get_top_level_ignored(self, psync: nf_core.pipelines.sync.PipelineSync) -> set[str]:
        """Helper function to get top-level part of relative directory of ignored files from psync.ignored_files."""
        top_level_ignored = {Path(f).parts[0] for f in psync.ignored_files}
        return top_level_ignored
