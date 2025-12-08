"""Tests for the download subcommand of nf-core tools"""

import os
import shutil
import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from unittest import mock

import pytest
import rich.progress_bar
import rich.table
import rich.text

from nf_core.pipelines.download import DownloadWorkflow
from nf_core.pipelines.download.docker import (
    DockerError,
    DockerFetcher,
    DockerProgress,
)

from ...utils import with_temporary_folder


#
# Test the DockerProgress subclass
#
class DockerProgressTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    def test_docker_progress_pull(self):
        with DockerProgress() as progress:
            assert progress.tasks == []
            progress.add_task(
                "Task 1", progress_type="docker", total=2, completed=1, current_log="example log", status="Pulling"
            )
            assert len(progress.tasks) == 1

            renderable = progress.get_renderable()
            assert isinstance(renderable, rich.console.Group), type(renderable)

            assert len(renderable.renderables) == 1
            table = renderable.renderables[0]
            assert isinstance(table, rich.table.Table)

            assert isinstance(table.columns[0]._cells[0], str)
            assert table.columns[0]._cells[0] == "[magenta]Task 1"
            assert isinstance(table.columns[2]._cells[0], str)
            assert table.columns[2]._cells[0] == "([blue]Pulling)"


#
# Test the DockerFetcher class
#
class DockerTest(unittest.TestCase):
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
    # Tests for 'pull_image'
    #
    # If Docker is installed, but the container can't be accessed because it does not exist or there are access
    # restrictions, a RuntimeWarning is raised due to the unavailability of the image.
    @pytest.mark.skipif(
        shutil.which("docker") is None,
        reason="Can't test what Docker does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.pipelines.download.docker.DockerProgress")
    @mock.patch("rich.progress.Task")
    def test_docker_pull_image_docker_installed(self, tmp_dir, mock_progress, mock_task):
        tmp_dir = Path(tmp_dir)
        docker_fetcher = DockerFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
        )
        docker_fetcher.progress = mock_progress()
        mock_task_obj = mock_task()

        # Test successful pull
        docker_fetcher.pull_image("hello-world", mock_task_obj)

        # Test successful pull with absolute URI (use tiny 3.5MB test container from the "Kogia" project: https://github.com/bschiffthaler/kogia)
        docker_fetcher.pull_image("docker.io/bschiffthaler/sed", mock_task_obj)

        # Test successful pull from wave
        docker_fetcher.pull_image(
            "community.wave.seqera.io/library/umi-transfer:1.0.0--e5b0c1a65b8173b6", mock_task_obj
        )

        # test image not found for several registries
        with pytest.raises(DockerError.ImageNotFoundError):
            docker_fetcher.pull_image("ghcr.io/not-a-real-registry/this-container-does-not-exist", mock_task_obj)

        with pytest.raises(DockerError.ImageNotFoundError):
            docker_fetcher.pull_image("docker.io/not-a-real-registry/this-container-does-not-exist", mock_task_obj)

        # test image not found for absolute URI.
        with pytest.raises(DockerError.ImageNotFoundError):
            docker_fetcher.pull_image("docker.io/bschiffthaler/nothingtopullhere", mock_task_obj)

        # Traffic from Github Actions to GitHub's Container Registry is unlimited, so no harm should be done here.
        with pytest.raises(DockerError.InvalidTagError):
            docker_fetcher.pull_image("ghcr.io/ewels/multiqc:go-rewrite", mock_task_obj)

    #
    # Tests for 'pull_and_save_image'
    #
    @pytest.mark.skipif(
        shutil.which("docker") is None,
        reason="Can't test what Docker does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.pipelines.download.docker.DockerProgress")
    def test_docker_pull_image_successfully(self, tmp_dir, mock_progress):
        tmp_dir = Path(tmp_dir)
        docker_fetcher = DockerFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
        )
        docker_fetcher.progress = mock_progress()
        docker_fetcher.pull_and_save_image("hello-world", tmp_dir / "hello-world.tar")

    #
    # Tests for 'save_image'
    #
    @pytest.mark.skipif(
        shutil.which("docker") is None,
        reason="Can't test what Docker does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.pipelines.download.docker.DockerProgress")
    @mock.patch("rich.progress.Task")
    def test_docker_save_image(self, tmp_dir, mock_progress, mock_task):
        tmp_dir = Path(tmp_dir)
        docker_fetcher = DockerFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
        )
        docker_fetcher.progress = mock_progress()
        mock_task_obj = mock_task()
        with pytest.raises(DockerError.ImageNotPulledError):
            docker_fetcher.save_image(
                "this-image-cannot-possibly-be-pulled-to-this-machine:latest",
                tmp_dir / "this-image-cannot-possibly-be-pulled-to-this-machine.tar",
                mock_task_obj,
            )

    #
    #
    # Tests for 'fetch_containers': this will test fetch remote containers automatically
    #
    @pytest.mark.skipif(
        shutil.which("singularity") is None and shutil.which("apptainer") is None,
        reason="Can't test what Singularity does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.utils.fetch_wf_config")
    def test_fetch_containers_docker(self, tmp_path, mock_fetch_wf_config):
        tmp_path = Path(tmp_path)
        download_obj = DownloadWorkflow(
            pipeline="dummy",
            outdir=tmp_path,
            container_library=None,
            container_system="docker",
        )
        download_obj.containers = [
            "helloworld",
            "helloooooooworld",
            "ewels/multiqc:gorewrite",
        ]
        # This list of fake container images should produce all kinds of ContainerErrors.
        # Test that they are all caught inside DockerFetcher.fetch_containers().
        docker_fetcher = DockerFetcher(
            outdir=tmp_path,
            container_library=download_obj.container_library,
            registry_set=download_obj.registry_set,
        )
        docker_fetcher.fetch_containers(
            download_obj.containers,
            download_obj.containers_remote,
            workflow_directory=Path("pipeline-dummy"),
        )

    #
    # Test for DockerFetcher.write_docker_load_command
    #
    @with_temporary_folder
    def test_docker_write_docker_load_message(self, tmp_path):
        tmp_path = Path(tmp_path)
        docker_fetcher = DockerFetcher(
            outdir=tmp_path,
            container_library=[],
            registry_set=[],
        )
        docker_output_dir = docker_fetcher.get_container_output_dir()
        docker_output_dir.mkdir()
        with redirect_stderr(StringIO()) as f:
            docker_fetcher.write_docker_load_message()

        # Check that the message looks ok
        assert str(docker_output_dir) in f.getvalue()
        assert "podman-load.sh" in f.getvalue()
        assert "docker-load.sh" in f.getvalue()

        # Check that the files were written and are executable
        assert (docker_output_dir / "docker-load.sh").exists()
        assert (docker_output_dir / "podman-load.sh").exists()
        assert os.access(docker_output_dir / "docker-load.sh", os.X_OK)
        assert os.access(docker_output_dir / "podman-load.sh", os.X_OK)
