"""Tests for the download subcommand of nf-core tools"""

import json
import logging
import os
import re
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pytest

import nf_core.pipelines.create.create
import nf_core.pipelines.download
import nf_core.pipelines.list
import nf_core.utils
from nf_core.pipelines.download import DownloadWorkflow
from nf_core.pipelines.download.workflow_repo import WorkflowRepo
from nf_core.synced_repo import SyncedRepo
from nf_core.utils import (
    NF_INSPECT_MIN_NF_VERSION,
    check_nextflow_version,
)

from ...utils import TEST_DATA_DIR, with_temporary_folder


class DownloadTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    @property
    def logged_levels(self) -> list[str]:
        return [record.levelname for record in self._caplog.records]

    @property
    def logged_messages(self) -> list[str]:
        return [record.message for record in self._caplog.records]

    def __contains__(self, item: str) -> bool:
        """Allows to check for log messages easily using the in operator inside a test:
        assert 'my log message' in self
        """
        return any(record.message == item for record in self._caplog.records if self._caplog)

    #
    # Tests for 'get_release_hash'
    #
    def test_get_release_hash_release_noauth(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        pipeline = "methylseq"

        try:
            # explicitly overwrite the gh_api authentication state
            # to force using the archive url
            from nf_core.utils import gh_api

            gh_api.lazy_init()
            gh_api.auth = None

            download_obj = DownloadWorkflow(pipeline=pipeline, revision="1.6")
            (
                download_obj.pipeline,
                download_obj.wf_revisions,
                download_obj.wf_branches,
            ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
            download_obj.get_revision_hash()
            assert download_obj.wf_sha[download_obj.revision[0]] == "b3e5e3b95aaf01d98391a62a10a3990c0a4de395"
            assert download_obj.outdir == Path("nf-core-methylseq_1.6")

            assert (
                download_obj.wf_download_url[download_obj.revision[0]]
                == "https://github.com/nf-core/methylseq/archive/b3e5e3b95aaf01d98391a62a10a3990c0a4de395.zip"
            )

        finally:
            gh_api.has_init = False
            gh_api.auth = None

    def test_get_release_hash_release(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        pipeline = "methylseq"
        download_obj = DownloadWorkflow(pipeline=pipeline, revision="1.6")
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
        download_obj.get_revision_hash()
        assert download_obj.wf_sha[download_obj.revision[0]] == "b3e5e3b95aaf01d98391a62a10a3990c0a4de395"
        assert download_obj.outdir == Path("nf-core-methylseq_1.6")
        assert (
            download_obj.wf_download_url[download_obj.revision[0]]
            == "https://api.github.com/repos/nf-core/methylseq/zipball/b3e5e3b95aaf01d98391a62a10a3990c0a4de395"
        )

    def test_get_release_hash_branch(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        # Exoseq pipeline is archived, so `dev` branch should be stable
        pipeline = "exoseq"
        download_obj = DownloadWorkflow(pipeline=pipeline, revision="dev")
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
        download_obj.get_revision_hash()
        assert download_obj.wf_sha[download_obj.revision[0]] == "819cbac792b76cf66c840b567ed0ee9a2f620db7"
        assert download_obj.outdir == Path("nf-core-exoseq_dev")
        assert (
            download_obj.wf_download_url[download_obj.revision[0]]
            == "https://api.github.com/repos/nf-core/exoseq/zipball/819cbac792b76cf66c840b567ed0ee9a2f620db7"
        )

    def test_get_release_hash_long_commit(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        # Exoseq pipeline is archived, so `dev` branch should be stable
        pipeline = "exoseq"
        revision = "819cbac792b76cf66c840b567ed0ee9a2f620db7"

        download_obj = DownloadWorkflow(pipeline=pipeline, revision=revision)
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
        download_obj.get_revision_hash()
        assert download_obj.wf_sha[download_obj.revision[0]] == revision
        assert download_obj.outdir == Path(f"nf-core-exoseq_{revision}")
        assert (
            download_obj.wf_download_url[download_obj.revision[0]]
            == f"https://api.github.com/repos/nf-core/exoseq/zipball/{revision}"
        )

    def test_get_release_hash_short_commit(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        # Exoseq pipeline is archived, so `dev` branch should be stable
        pipeline = "exoseq"
        revision = "819cbac792b76cf66c840b567ed0ee9a2f620db7"
        short_rev = revision[:7]

        download_obj = DownloadWorkflow(pipeline="exoseq", revision=short_rev)
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
        download_obj.get_revision_hash()
        print(download_obj)
        assert download_obj.wf_sha[download_obj.revision[0]] == revision
        assert download_obj.outdir == Path(f"nf-core-exoseq_{short_rev}")
        assert (
            download_obj.wf_download_url[download_obj.revision[0]]
            == f"https://api.github.com/repos/nf-core/exoseq/zipball/{revision}"
        )

    def test_get_release_hash_non_existent_release(self):
        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        pipeline = "methylseq"
        download_obj = DownloadWorkflow(pipeline=pipeline, revision="thisisfake")
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(pipeline, wfs)
        with pytest.raises(AssertionError):
            download_obj.get_revision_hash()

    #
    # Tests for 'download_wf_files'
    #
    @with_temporary_folder
    def test_download_wf_files(self, outdir):
        outdir = Path(outdir)
        download_obj = DownloadWorkflow(pipeline="nf-core/methylseq", revision="1.6")
        download_obj.outdir = outdir
        download_obj.wf_sha = {"1.6": "b3e5e3b95aaf01d98391a62a10a3990c0a4de395"}
        download_obj.wf_download_url = {
            "1.6": "https://api.github.com/repos/nf-core/methylseq/zipball/b3e5e3b95aaf01d98391a62a10a3990c0a4de395"
        }
        rev = download_obj.download_wf_files(
            download_obj.revision[0],
            download_obj.wf_sha[download_obj.revision[0]],
            download_obj.wf_download_url[download_obj.revision[0]],
        )

        assert ((outdir / rev) / "main.nf").exists()

    #
    # Tests for 'download_configs'
    #
    @with_temporary_folder
    def test_download_configs(self, outdir):
        outdir = Path(outdir)
        download_obj = DownloadWorkflow(pipeline="nf-core/methylseq", revision="1.6")
        download_obj.outdir = outdir
        download_obj.download_configs()
        assert (outdir / "configs") / "nfcore_custom.config"

    #
    # Tests for 'wf_use_local_configs'
    #
    @with_temporary_folder
    def test_wf_use_local_configs(self, tmp_path):
        tmp_path = Path(tmp_path)
        # Get a workflow and configs
        test_pipeline_dir = tmp_path / "nf-core-testpipeline"
        create_obj = nf_core.pipelines.create.create.PipelineCreate(
            "testpipeline",
            "This is a test pipeline",
            "Test McTestFace",
            no_git=True,
            outdir=test_pipeline_dir,
        )
        create_obj.init_pipeline()

        with tempfile.TemporaryDirectory() as test_outdir:
            download_obj = DownloadWorkflow(pipeline="dummy", revision="1.2.0", outdir=test_outdir)
            shutil.copytree(test_pipeline_dir, Path(test_outdir, "workflow"))
            download_obj.download_configs()

            # Test the function
            download_obj.wf_use_local_configs("workflow")
            wf_config = nf_core.utils.fetch_wf_config(Path(test_outdir, "workflow"), cache_config=False)
            assert wf_config["params.custom_config_base"] == f"{test_outdir}/workflow/../configs/"

    #
    # Test that `find_container_images` (uses `nextflow inspect`) fetches the correct Docker images
    #
    @pytest.mark.skipif(
        shutil.which("nextflow") is None or not check_nextflow_version(NF_INSPECT_MIN_NF_VERSION),
        reason="Can't run test that requires nextflow to run if not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.utils.fetch_wf_config")
    def test_containers_pipeline_singularity(self, tmp_path, mock_fetch_wf_config):
        tmp_path = Path(tmp_path)
        assert check_nextflow_version(NF_INSPECT_MIN_NF_VERSION) is True

        # Set up test
        container_system = "singularity"
        mock_pipeline_dir = TEST_DATA_DIR / "mock_pipeline_containers"
        refererence_json_dir = mock_pipeline_dir / "per_profile_output"
        # First check that `-profile singularity` produces the same output as the reference
        download_obj = DownloadWorkflow(pipeline="dummy", outdir=tmp_path, container_system=container_system)
        mock_fetch_wf_config.return_value = {}

        # Run get containers with `nextflow inspect`
        entrypoint = "main_passing_test.nf"
        download_obj.find_container_images(mock_pipeline_dir, "dummy-revision", entrypoint=entrypoint)

        # Store the containers found by the new method
        found_containers = set(download_obj.containers)

        # Load the reference containers
        with open(refererence_json_dir / f"{container_system}_containers.json") as fh:
            ref_containers = json.load(fh)
            ref_container_strs = set(ref_containers.values())

        # Now check that they contain the same containers
        assert found_containers == ref_container_strs, (
            f"Containers found in pipeline by `nextflow inspect`: {found_containers}\n"
            f"Containers that should've been found: {ref_container_strs}"
        )

    #
    # Test that `find_container_images` (uses `nextflow inspect`) fetches the correct Docker images
    #
    @pytest.mark.skipif(
        shutil.which("nextflow") is None or not check_nextflow_version(NF_INSPECT_MIN_NF_VERSION),
        reason=f"Can't run test that requires Nextflow >= {NF_INSPECT_MIN_NF_VERSION} to run if not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.utils.fetch_wf_config")
    def test_containers_pipeline_docker(self, tmp_path, mock_fetch_wf_config):
        tmp_path = Path(tmp_path)
        assert check_nextflow_version(NF_INSPECT_MIN_NF_VERSION) is True

        # Set up test
        container_system = "docker"
        mock_pipeline_dir = TEST_DATA_DIR / "mock_pipeline_containers"
        refererence_json_dir = mock_pipeline_dir / "per_profile_output"
        # First check that `-profile singularity` produces the same output as the reference
        download_obj = DownloadWorkflow(pipeline="dummy", outdir=tmp_path, container_system=container_system)
        mock_fetch_wf_config.return_value = {}

        # Run get containers with `nextflow inspect`
        entrypoint = "main_passing_test.nf"
        download_obj.find_container_images(mock_pipeline_dir, "dummy-revision", entrypoint=entrypoint)

        # Store the containers found by the new method
        found_containers = set(download_obj.containers)

        # Load the reference containers
        with open(refererence_json_dir / f"{container_system}_containers.json") as fh:
            ref_containers = json.load(fh)
            ref_container_strs = set(ref_containers.values())

        # Now check that they contain the same containers
        assert found_containers == ref_container_strs, (
            f"Containers found in pipeline by `nextflow inspect`: {found_containers}\n"
            f"Containers that should've been found: {ref_container_strs}"
        )

    #
    # Tests for the main entry method 'download_workflow'
    #

    # We do not want to download all containers, so we mock the download by just touching the singularity files
    def mock_download_file(self, remote_path: str, output_path: str):
        Path(output_path).touch()  # Create an empty file at the output path

    @with_temporary_folder
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.check_and_set_implementation"
    )  # This is to make sure that we do not check for Singularity/Apptainer installation
    @mock.patch.object(nf_core.pipelines.download.singularity.FileDownloader, "download_file", new=mock_download_file)
    def test_download_workflow_with_success(self, tmp_dir, mock_check_and_set_implementation):
        tmp_dir = Path(tmp_dir)
        os.environ["NXF_SINGULARITY_CACHEDIR"] = str(tmp_dir / "foo")

        download_obj = DownloadWorkflow(
            pipeline="nf-core/bamtofastq",
            outdir=tmp_dir / "new",
            container_system="singularity",
            revision="2.2.0",
            compress_type="none",
            container_cache_utilisation="copy",
            parallel=1,
        )

        download_obj.include_configs = True  # suppress prompt, because stderr.is_interactive doesn't.
        download_obj.download_workflow()

    #
    # Test Download for Seqera Platform
    #
    @with_temporary_folder
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.check_and_set_implementation"
    )  # This is to make sure that we do not check for Singularity/Apptainer installation
    @mock.patch("nf_core.pipelines.download.singularity.SingularityFetcher.fetch_containers")
    def test_download_workflow_for_platform(
        self,
        tmp_dir,
        mock_fetch_containers,
        mock_check_and_set_implementation,
    ):
        tmp_dir = Path(tmp_dir)
        download_obj = DownloadWorkflow(
            pipeline="nf-core/rnaseq",
            revision=("3.19.0", "3.17.0"),
            compress_type="none",
            platform=True,
            container_system="singularity",
        )

        download_obj.include_configs = False  # suppress prompt, because stderr.is_interactive doesn't.

        assert isinstance(download_obj.revision, list) and len(download_obj.revision) == 2
        assert isinstance(download_obj.wf_sha, dict) and len(download_obj.wf_sha) == 0
        assert isinstance(download_obj.wf_download_url, dict) and len(download_obj.wf_download_url) == 0

        wfs = nf_core.pipelines.list.Workflows()
        wfs.get_remote_workflows()
        (
            download_obj.pipeline,
            download_obj.wf_revisions,
            download_obj.wf_branches,
        ) = nf_core.utils.get_repo_releases_branches(download_obj.pipeline, wfs)

        download_obj.get_revision_hash()

        # download_obj.wf_download_url is not set for Seqera Platform downloads, but the sha values are
        assert isinstance(download_obj.wf_sha, dict) and len(download_obj.wf_sha) == 2
        assert isinstance(download_obj.wf_download_url, dict) and len(download_obj.wf_download_url) == 0

        # The outdir for multiple revisions is the pipeline name and date: e.g. nf-core-rnaseq_2023-04-27_18-54
        assert isinstance(download_obj.outdir, Path)
        assert bool(re.search(r"nf-core-rnaseq_\d{4}-\d{2}-\d{1,2}_\d{1,2}-\d{1,2}", str(download_obj.outdir), re.S))

        download_obj.output_filename = download_obj.outdir.with_suffix(".git")
        download_obj.download_workflow_platform(location=tmp_dir)

        assert download_obj.workflow_repo
        assert isinstance(download_obj.workflow_repo, WorkflowRepo)
        assert issubclass(type(download_obj.workflow_repo), SyncedRepo)

        # corroborate that the other revisions are inaccessible to the user.
        all_tags = {tag.name for tag in download_obj.workflow_repo.tags}
        all_heads = {head.name for head in download_obj.workflow_repo.heads}

        assert set(download_obj.revision) == all_tags
        # assert that the download has a "latest" branch.
        assert "latest" in all_heads

        # download_obj.download_workflow_platform(location=tmp_dir) will run `nextflow inspect` for each revision
        # This means that the containers in download_obj.containers are the containers the last specified revision i.e. 3.17
        assert isinstance(download_obj.containers, list) and len(download_obj.containers) == 39
        assert (
            "https://depot.galaxyproject.org/singularity/bbmap:39.10--h92535d8_0" in download_obj.containers
        )  # direct definition

        # clean-up
        # remove "nf-core-rnaseq*" directories
        for path in Path().cwd().glob("nf-core-rnaseq*"):
            shutil.rmtree(path)

    #
    # Brief test adding a single custom tag to Seqera Platform download
    #
    @mock.patch("nf_core.pipelines.download.singularity.SingularityFetcher.fetch_containers")
    @with_temporary_folder
    def test_download_workflow_for_platform_with_one_custom_tag(self, _, tmp_dir):
        tmp_dir = Path(tmp_dir)
        download_obj = DownloadWorkflow(
            pipeline="nf-core/rnaseq",
            revision=("3.9"),
            compress_type="none",
            platform=True,
            container_system=None,
            additional_tags=("3.9=cool_revision",),
        )
        assert isinstance(download_obj.additional_tags, list) and len(download_obj.additional_tags) == 1

        # clean-up
        # remove "nf-core-rnaseq*" directories
        for path in Path().cwd().glob("nf-core-rnaseq*"):
            shutil.rmtree(path)

    #
    # Test adding custom tags to Seqera Platform download (full test)
    #
    @mock.patch("nf_core.pipelines.download.singularity.SingularityFetcher.fetch_containers")
    @with_temporary_folder
    def test_download_workflow_for_platform_with_custom_tags(self, _, tmp_dir):
        tmp_dir = Path(tmp_dir)
        with self._caplog.at_level(logging.INFO):
            from git.refs.tag import TagReference

            download_obj = DownloadWorkflow(
                pipeline="nf-core/rnaseq",
                revision=("3.7", "3.9"),
                compress_type="none",
                platform=True,
                container_system=None,
                additional_tags=(
                    "3.7=a.tad.outdated",
                    "3.9=cool_revision",
                    "3.9=invalid tag",
                    "3.14.0=not_included",
                    "What is this?",
                ),
            )

            download_obj.include_configs = False  # suppress prompt, because stderr.is_interactive doesn't.

            assert isinstance(download_obj.revision, list) and len(download_obj.revision) == 2
            assert isinstance(download_obj.wf_sha, dict) and len(download_obj.wf_sha) == 0
            assert isinstance(download_obj.wf_download_url, dict) and len(download_obj.wf_download_url) == 0
            assert isinstance(download_obj.additional_tags, list) and len(download_obj.additional_tags) == 5

            wfs = nf_core.pipelines.list.Workflows()
            wfs.get_remote_workflows()
            (
                download_obj.pipeline,
                download_obj.wf_revisions,
                download_obj.wf_branches,
            ) = nf_core.utils.get_repo_releases_branches(download_obj.pipeline, wfs)

            download_obj.get_revision_hash()
            download_obj.output_filename = f"{download_obj.outdir}.git"
            download_obj.download_workflow_platform(location=tmp_dir)

            assert download_obj.workflow_repo
            assert isinstance(download_obj.workflow_repo, WorkflowRepo)
            assert issubclass(type(download_obj.workflow_repo), SyncedRepo)
            assert "Locally cached repository: nf-core/rnaseq, revisions 3.7, 3.9" in repr(download_obj.workflow_repo)

            # assert that every additional tag has been passed on to the WorkflowRepo instance
            assert download_obj.additional_tags == download_obj.workflow_repo.additional_tags

            # assert that the additional tags are all TagReference objects
            assert all(isinstance(tag, TagReference) for tag in download_obj.workflow_repo.tags)

            workflow_repo_tags = {tag.name for tag in download_obj.workflow_repo.tags}
            assert len(workflow_repo_tags) == 4
            # the invalid/malformed additional_tags should not have been added.
            assert all(tag in workflow_repo_tags for tag in {"3.7", "a.tad.outdated", "cool_revision", "3.9"})
            assert not any(tag in workflow_repo_tags for tag in {"invalid tag", "not_included", "What is this?"})

            assert all(
                log in self.logged_messages
                for log in {
                    "[red]Could not apply invalid `--tag` specification[/]: '3.9=invalid tag'",
                    "[red]Adding tag 'not_included' to '3.14.0' failed.[/]\n Mind that '3.14.0' must be a valid git reference that resolves to a commit.",
                    "[red]Could not apply invalid `--tag` specification[/]: 'What is this?'",
                }
            )

            # clean-up
            # remove "nf-core-rnaseq*" directories
            for path in Path().cwd().glob("nf-core-rnaseq*"):
                shutil.rmtree(path)
