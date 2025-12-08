"""Tests for the download subcommand of nf-core tools"""

import logging
import os
import shutil
import unittest
from pathlib import Path
from unittest import mock

import pytest
import requests
import rich.progress_bar
import rich.table
import rich.text

from nf_core.pipelines.download import DownloadWorkflow
from nf_core.pipelines.download.container_fetcher import ContainerProgress
from nf_core.pipelines.download.singularity import (
    FileDownloader,
    SingularityError,
    SingularityFetcher,
    SingularityProgress,
)
from nf_core.pipelines.download.utils import DownloadError

from ...utils import TEST_DATA_DIR, with_temporary_folder


#
# Test the SingularityProgress subclass
#
class SingularityProgressTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    # Test the "singularity_pull" progress type
    def test_singularity_pull_progress(self):
        with SingularityProgress() as progress:
            assert progress.tasks == []
            progress.add_task(
                "Task 1", progress_type="singularity_pull", total=42, completed=11, current_log="example log"
            )
            assert len(progress.tasks) == 1

            renderable = progress.get_renderable()
            assert isinstance(renderable, rich.console.Group), type(renderable)

            assert len(renderable.renderables) == 1
            table = renderable.renderables[0]
            assert isinstance(table, rich.table.Table)

            assert isinstance(table.columns[0]._cells[0], str)
            assert table.columns[0]._cells[0] == "[magenta]Task 1"

            assert isinstance(table.columns[1]._cells[0], str)
            assert table.columns[1]._cells[0] == "[blue]example log"

            assert isinstance(table.columns[2]._cells[0], rich.progress_bar.ProgressBar)
            assert table.columns[2]._cells[0].completed == 11
            assert table.columns[2]._cells[0].total == 42

    # Test the "download" progress type
    def test_singularity_download_progress(self):
        with SingularityProgress() as progress:
            assert progress.tasks == []
            progress.add_task("Task 1", progress_type="download", total=42, completed=11)
            assert len(progress.tasks) == 1

            renderable = progress.get_renderable()
            assert isinstance(renderable, rich.console.Group), type(renderable)

            assert len(renderable.renderables) == 1
            table = renderable.renderables[0]
            assert isinstance(table, rich.table.Table)

            assert isinstance(table.columns[0]._cells[0], str)
            assert table.columns[0]._cells[0] == "[blue]Task 1"

            assert isinstance(table.columns[1]._cells[0], rich.progress_bar.ProgressBar)
            assert table.columns[1]._cells[0].completed == 11
            assert table.columns[1]._cells[0].total == 42

            assert isinstance(table.columns[2]._cells[0], str)
            assert table.columns[2]._cells[0] == "[progress.percentage]26.2%"

            assert isinstance(table.columns[3]._cells[0], str)
            assert table.columns[3]._cells[0] == "•"

            assert isinstance(table.columns[4]._cells[0], rich.text.Text)
            assert table.columns[4]._cells[0]._text == ["11/42 bytes"]

            assert isinstance(table.columns[5]._cells[0], str)
            assert table.columns[5]._cells[0] == "•"

            assert isinstance(table.columns[6]._cells[0], rich.text.Text)
            assert table.columns[6]._cells[0]._text == ["?"]


#
# Test the SingularityFetcher class
#
class SingularityTest(unittest.TestCase):
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
    # Tests for 'singularity_pull_image'
    #
    # If Singularity is installed, but the container can't be accessed because it does not exist or there are access
    # restrictions, a RuntimeWarning is raised due to the unavailability of the image.
    @pytest.mark.skipif(
        shutil.which("singularity") is None and shutil.which("apptainer") is None,
        reason="Can't test what Singularity does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.pipelines.download.singularity.SingularityProgress")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_singularity_pull_image_singularity_installed(self, tmp_dir, mock_cachedir_prompt, mock_progress):
        tmp_dir = Path(tmp_dir)
        singularity_fetcher = SingularityFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
            container_cache_utilisation="none",
            container_cache_index=None,
        )
        singularity_fetcher.check_and_set_implementation()
        singularity_fetcher.progress = mock_progress()
        singularity_fetcher.registry_set = {}
        # Test successful pull
        assert singularity_fetcher.pull_image("hello-world", tmp_dir / "hello-world.sif", "docker.io") is True

        # Pull again, but now the image already exists
        assert singularity_fetcher.pull_image("hello-world", tmp_dir / "hello-world.sif", "docker.io") is False

        # Test successful pull with absolute URI (use tiny 3.5MB test container from the "Kogia" project: https://github.com/bschiffthaler/kogia)
        assert singularity_fetcher.pull_image("docker.io/bschiffthaler/sed", tmp_dir / "sed.sif", "docker.io") is True

        # Test successful pull with absolute oras:// URI
        assert (
            singularity_fetcher.pull_image(
                "oras://community.wave.seqera.io/library/umi-transfer:1.0.0--e5b0c1a65b8173b6",
                tmp_dir / "umi-transfer-oras.sif",
                "docker.io",
            )
            is True
        )

        # try pulling Docker container image with oras://
        with pytest.raises(SingularityError.NoSingularityContainerError):
            singularity_fetcher.pull_image(
                "oras://ghcr.io/matthiaszepper/umi-transfer:dev",
                tmp_dir / "umi-transfer-oras_impostor.sif",
                "docker.io",
            )

        # try to pull from non-existing registry (Name change hello-world_new.sif is needed, otherwise ImageExistsError is raised before attempting to pull.)
        with pytest.raises(SingularityError.RegistryNotFoundError):
            singularity_fetcher.pull_image(
                "hello-world",
                tmp_dir / "break_the_registry_test.sif",
                "register-this-domain-to-break-the-test.io",
            )

        # test Image not found for several registries
        with pytest.raises(SingularityError.ImageNotFoundError):
            singularity_fetcher.pull_image("a-container", tmp_dir / "acontainer.sif", "quay.io")

        with pytest.raises(SingularityError.ImageNotFoundError):
            singularity_fetcher.pull_image("a-container", tmp_dir / "acontainer.sif", "docker.io")

        with pytest.raises(SingularityError.ImageNotFoundError):
            singularity_fetcher.pull_image("a-container", tmp_dir / "acontainer.sif", "ghcr.io")

        # test Image not found for absolute URI.
        with pytest.raises(SingularityError.ImageNotFoundError):
            singularity_fetcher.pull_image(
                "docker.io/bschiffthaler/nothingtopullhere",
                tmp_dir / "nothingtopullhere.sif",
                "docker.io",
            )

        # Traffic from Github Actions to GitHub's Container Registry is unlimited, so no harm should be done here.
        with pytest.raises(SingularityError.InvalidTagError):
            singularity_fetcher.pull_image(
                "ewels/multiqc:go-rewrite",
                tmp_dir / "multiqc-go.sif",
                "ghcr.io",
            )

    #
    # Tests for 'fetch_containers': this will test fetch remote containers automatically
    #
    @pytest.mark.skipif(
        shutil.which("singularity") is None and shutil.which("apptainer") is None,
        reason="Can't test what Singularity does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.utils.fetch_wf_config")
    @mock.patch("nf_core.pipelines.download.singularity.SingularityFetcher.gather_registries")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_fetch_containers_singularity(
        self, tmp_path, mock_cachedir_prompt, mock_gather_registries, mock_fetch_wf_config
    ):
        tmp_path = Path(tmp_path)
        download_obj = DownloadWorkflow(
            pipeline="dummy",
            outdir=tmp_path,
            container_library=("mirage-the-imaginative-registry.io", "quay.io", "ghcr.io", "docker.io"),
            container_system="singularity",
        )
        download_obj.containers = [
            "helloworld",
            "helloooooooworld",
            "ewels/multiqc:gorewrite",
        ]
        assert len(download_obj.container_library) == 4
        # This list of fake container images should produce all kinds of ContainerErrors.
        # Test that they are all caught inside SingularityFetcher.fetch_containers().
        singularity_fetcher = SingularityFetcher(
            outdir=tmp_path,
            container_library=download_obj.container_library,
            registry_set=download_obj.registry_set,
            container_cache_utilisation="none",
            container_cache_index=None,
        )
        singularity_fetcher.fetch_containers(
            download_obj.containers,
            download_obj.containers_remote,
            workflow_directory=Path("pipeline-dummy"),
        )

    #
    # Tests for the 'symlink_registries' function
    #

    # Simple file name with no registry in it
    @with_temporary_folder
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.check_and_set_implementation"
    )  # This is to make sure that we do not check for Singularity/Apptainer installation
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    @mock.patch("pathlib.Path.mkdir")
    @mock.patch("pathlib.Path.symlink_to")
    @mock.patch("os.symlink")
    @mock.patch("os.open")
    @mock.patch("os.close")
    @mock.patch("pathlib.Path.name")
    @mock.patch("pathlib.Path.parent")
    def test_symlink_singularity_images(
        self,
        tmp_path,
        mock_dirname,
        mock_basename,
        mock_close,
        mock_open,
        mock_os_symlink,
        mock_symlink,
        mock_makedirs,
        mock_prompt_singularity_cachedir_creation,
        mock_check_and_set_implementation,
    ):
        # Setup
        tmp_path = Path(tmp_path)
        with (
            mock.patch.object(Path, "name", new_callable=mock.PropertyMock) as mock_basename,
            mock.patch.object(Path, "parent", new_callable=mock.PropertyMock) as mock_dirname,
        ):
            mock_dirname.return_value = tmp_path / "path/to"
            mock_basename.return_value = "singularity-image.img"
            mock_open.return_value = 12  # file descriptor
            mock_close.return_value = 12  # file descriptor
            mock_prompt_singularity_cachedir_creation.return_value = False

            registries = [
                "quay.io",
                "community-cr-prod.seqera.io/docker/registry/v2",
                "depot.galaxyproject.org/singularity",
            ]
            fetcher = SingularityFetcher(
                outdir=tmp_path,
                container_library=[],
                registry_set=registries,
                container_cache_utilisation="none",
                container_cache_index=None,
            )
            fetcher.registry_set = registries
            fetcher.symlink_registries(tmp_path / "path/to/singularity-image.img")

            # Check that os.makedirs was called with the correct arguments
            mock_makedirs.assert_any_call(exist_ok=True)

            # Check that os.open was called with the correct arguments
            mock_open.assert_any_call(tmp_path / "path/to", os.O_RDONLY)

            # Check that os.symlink was called with the correct arguments
            expected_calls = [
                mock.call(
                    Path("./singularity-image.img"),
                    Path("./quay.io-singularity-image.img"),
                    dir_fd=12,
                ),
                mock.call(
                    Path("./singularity-image.img"),
                    Path("./community-cr-prod.seqera.io-docker-registry-v2-singularity-image.img"),
                    dir_fd=12,
                ),
                mock.call(
                    Path("./singularity-image.img"),
                    Path("./depot.galaxyproject.org-singularity-singularity-image.img"),
                    dir_fd=12,
                ),
            ]
            mock_os_symlink.assert_has_calls(expected_calls, any_order=True)

    # File name with registry in it
    @with_temporary_folder
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.check_and_set_implementation"
    )  # This is to make sure that we do not check for Singularity/Apptainer installation
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    @mock.patch("pathlib.Path.mkdir")
    @mock.patch("pathlib.Path.symlink_to")
    @mock.patch("os.symlink")
    @mock.patch("os.open")
    @mock.patch("os.close")
    @mock.patch("re.sub")
    @mock.patch("pathlib.Path.name")
    @mock.patch("pathlib.Path.parent")
    def test_symlink_singularity_symlink_registries(
        self,
        tmp_path,
        mock_dirname,
        mock_basename,
        mock_resub,
        mock_close,
        mock_open,
        mock_os_symlink,
        mock_symlink,
        mock_makedirs,
        mock_prompt_singularity_cachedir_creation,
        mock_check_and_set_implementation,
    ):
        tmp_path = Path(tmp_path)
        # Setup
        with (
            mock.patch.object(Path, "name", new_callable=mock.PropertyMock) as mock_basename,
            mock.patch.object(Path, "parent", new_callable=mock.PropertyMock) as mock_dirname,
        ):
            mock_resub.return_value = "singularity-image.img"
            mock_dirname.return_value = tmp_path / "path/to"
            mock_basename.return_value = "quay.io-singularity-image.img"
            mock_open.return_value = 12  # file descriptor
            mock_close.return_value = 12  # file descriptor
            mock_prompt_singularity_cachedir_creation.return_value = False

            # Call the method with registry name included - should not happen, but preserve it then.

            registries = [
                "quay.io",  # Same as in the filename
                "community-cr-prod.seqera.io/docker/registry/v2",
            ]
            fetcher = SingularityFetcher(
                outdir=tmp_path,
                container_library=[],
                registry_set=registries,
                container_cache_utilisation="none",
                container_cache_index=None,
            )
            fetcher.registry_set = registries
            fetcher.symlink_registries(tmp_path / "path/to/quay.io-singularity-image.img")

            # Check that os.makedirs was called with the correct arguments
            mock_makedirs.assert_called_once_with(exist_ok=True)

            # Check that os.symlink was called with the correct arguments
            # assert_called_once_with also tells us that there was no attempt to
            # - symlink to itself
            # - symlink to the same registry
            mock_os_symlink.assert_called_once_with(
                Path("./quay.io-singularity-image.img"),
                Path(
                    "./community-cr-prod.seqera.io-docker-registry-v2-singularity-image.img"
                ),  # "quay.io-" has been trimmed
                dir_fd=12,
            )

            # Normally it would be called for each registry, but since quay.io is part of the name, it
            # will only be called once, as no symlink to itself must be created.
            mock_open.assert_called_once_with(tmp_path / "path/to", os.O_RDONLY)

    #
    # Tests for 'pull_image'
    #
    @pytest.mark.skipif(
        shutil.which("singularity") is None and shutil.which("apptainer") is None,
        reason="Can't test what Singularity does if it's not installed.",
    )
    @with_temporary_folder
    @mock.patch("nf_core.pipelines.download.singularity.SingularityProgress")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_singularity_pull_image_successfully(self, tmp_dir, mock_cachedir_prompt, mock_progress):
        tmp_dir = Path(tmp_dir)
        singularity_fetcher = SingularityFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
            container_cache_utilisation="none",
            container_cache_index=None,
        )
        singularity_fetcher.registry_set = set()
        singularity_fetcher.check_and_set_implementation()
        singularity_fetcher.progress = mock_progress()
        singularity_fetcher.pull_image("hello-world", tmp_dir / "yet-another-hello-world.sif", "docker.io")

    #
    @with_temporary_folder
    @mock.patch("nf_core.utils.fetch_wf_config")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_gather_registries_singularity(self, tmp_path, mock_cachedir_prompt, mock_fetch_wf_config):
        tmp_path = Path(tmp_path)
        container_library = ["quay.io"]
        singularity_fetcher = SingularityFetcher(
            outdir=tmp_path,
            container_library=container_library,
            registry_set=container_library,
            container_cache_utilisation="none",
            container_cache_index=None,
        )
        mock_fetch_wf_config.return_value = {
            "apptainer.registry": "apptainer-registry.io",
            "docker.registry": "docker.io",
            "podman.registry": "podman-registry.io",
            "singularity.registry": "singularity-registry.io",
            "someother.registry": "fake-registry.io",
        }
        singularity_fetcher.registry_set = singularity_fetcher.gather_registries(tmp_path)
        assert singularity_fetcher.registry_set
        assert isinstance(singularity_fetcher.registry_set, set)
        assert len(singularity_fetcher.registry_set) == 8

        assert "quay.io" in singularity_fetcher.registry_set  # default registry, if no container library is provided.
        assert (
            "depot.galaxyproject.org/singularity" in singularity_fetcher.registry_set
        )  # default registry, often hardcoded in modules
        assert "community.wave.seqera.io/library" in singularity_fetcher.registry_set  # Seqera containers Docker
        assert (
            "community-cr-prod.seqera.io/docker/registry/v2" in singularity_fetcher.registry_set
        )  # Seqera containers Singularity https:// download
        assert "apptainer-registry.io" in singularity_fetcher.registry_set
        assert "docker.io" in singularity_fetcher.registry_set
        assert "podman-registry.io" in singularity_fetcher.registry_set
        assert "singularity-registry.io" in singularity_fetcher.registry_set
        # it should only pull the apptainer, docker, podman and singularity registry from the config, but not any registry.
        assert "fake-registry.io" not in singularity_fetcher.registry_set

    #
    # If Singularity is not installed, it raises a OSError because the singularity command can't be found.
    #
    @pytest.mark.skipif(
        shutil.which("singularity") is not None or shutil.which("apptainer") is not None,
        reason="Can't test how the code behaves when singularity is not installed if it is.",
    )
    @with_temporary_folder
    @mock.patch("rich.progress.Progress.add_task")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_singularity_pull_image_singularity_not_installed(self, tmp_dir, mock_rich_progress, mock_cachedir_prompt):
        tmp_dir = Path(tmp_dir)
        fetcher = SingularityFetcher(
            outdir=tmp_dir,
            container_library=[],
            registry_set=[],
            container_cache_utilisation="none",
            container_cache_index=None,
        )
        with pytest.raises(OSError):
            fetcher.check_and_set_implementation()

    #
    # Test for 'get_container_filename' function
    #

    @mock.patch("nf_core.pipelines.download.singularity.SingularityFetcher.check_and_set_implementation")
    @mock.patch(
        "nf_core.pipelines.download.singularity.SingularityFetcher.prompt_singularity_cachedir_creation"
    )  # This is to make sure that we do not prompt for a Singularity cachedir
    def test_singularity_get_container_filename(self, mock_cachedir_prompt, mock_check_and_set_implementation):
        registries = [
            "docker.io",
            "quay.io",
            "depot.galaxyproject.org/singularity",
            "community.wave.seqera.io/library",
            "community-cr-prod.seqera.io/docker/registry/v2",
        ]

        fetcher = SingularityFetcher(
            outdir=Path("test_singularity_get_container_filename"),
            container_library=[],
            registry_set=registries,
            container_cache_utilisation="none",
            container_cache_index=None,
        )

        fetcher.registry_set = registries
        # Test --- galaxy URL #
        result = fetcher.get_container_filename(
            "https://depot.galaxyproject.org/singularity/bbmap:38.93--he522d1c_0",
        )
        assert result == "bbmap-38.93--he522d1c_0.img"

        # Test --- mulled containers #
        result = fetcher.get_container_filename(
            "quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:59cdd445419f14abac76b31dd0d71217994cbcc9-0",
        )
        assert (
            result
            == "biocontainers-mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2-59cdd445419f14abac76b31dd0d71217994cbcc9-0.img"
        )

        # Test --- Docker containers without registry #
        result = fetcher.get_container_filename("nf-core/ubuntu:20.04")
        assert result == "nf-core-ubuntu-20.04.img"

        # Test --- Docker container with explicit registry -> should be trimmed #
        result = fetcher.get_container_filename("docker.io/nf-core/ubuntu:20.04")
        assert result == "nf-core-ubuntu-20.04.img"

        # Test --- Docker container with explicit registry not in registry list -> can't be trimmed
        result = fetcher.get_container_filename("mirage-the-imaginative-registry.io/nf-core/ubuntu:20.04")
        assert result == "mirage-the-imaginative-registry.io-nf-core-ubuntu-20.04.img"

        # Test --- Seqera Docker containers: Trimmed, because it is hard-coded in the registry set.
        result = fetcher.get_container_filename("community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264")
        assert result == "coreutils-9.5--ae99c88a9b28c264.img"

        # Test --- Seqera Singularity containers: Trimmed, because it is hard-coded in the registry set.
        result = fetcher.get_container_filename(
            "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data",
        )
        assert result == "blobs-sha256-c2-c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975-data.img"

        # Test --- Seqera Oras containers: Trimmed, because it is hard-coded in the registry set.
        result = fetcher.get_container_filename(
            "oras://community.wave.seqera.io/library/umi-transfer:1.0.0--e5b0c1a65b8173b6",
        )
        assert result == "umi-transfer-1.0.0--e5b0c1a65b8173b6.img"

        # Test --- SIF Singularity container with explicit registry -> should be trimmed #
        result = fetcher.get_container_filename(
            "docker.io-hashicorp-vault-1.16-sha256:e139ff28c23e1f22a6e325696318141259b177097d8e238a3a4c5b84862fadd8.sif",
        )
        assert (
            result == "hashicorp-vault-1.16-sha256-e139ff28c23e1f22a6e325696318141259b177097d8e238a3a4c5b84862fadd8.sif"
        )

        # Test --- SIF Singularity container without registry #
        result = fetcher.get_container_filename(
            "singularity-hpc/shpc/tests/testdata/salad_latest.sif",
        )
        assert result == "singularity-hpc-shpc-tests-testdata-salad_latest.sif"

        # Test --- Singularity container from a Singularity registry (and version tag) #
        result = fetcher.get_container_filename(
            "library://pditommaso/foo/bar.sif:latest",
        )
        assert result == "pditommaso-foo-bar-latest.sif"

        # Test --- galaxy URL but no registry given #
        fetcher.registry_set = []
        result = fetcher.get_container_filename("https://depot.galaxyproject.org/singularity/bbmap:38.93--he522d1c_0")
        assert result == "depot.galaxyproject.org-singularity-bbmap-38.93--he522d1c_0.img"

    #
    # Test for '--singularity-cache remote --singularity-cache-index'. Provide a list of containers already available in a remote location.
    #
    @with_temporary_folder
    def test_remote_container_functionality(self, tmp_dir):
        tmp_dir = Path(tmp_dir)
        os.environ["NXF_SINGULARITY_CACHEDIR"] = str(tmp_dir / "foo")

        download_obj = DownloadWorkflow(
            pipeline="nf-core/rnaseq",
            outdir=(tmp_dir / "new"),
            revision="3.9",
            compress_type="none",
            container_cache_index=Path(TEST_DATA_DIR, "testdata_remote_containers.txt"),
            container_system="singularity",
        )

        download_obj.include_configs = False  # suppress prompt, because stderr.is_interactive doesn't.

        # test if the settings are changed to mandatory defaults, if an external cache index is used.
        assert download_obj.container_cache_utilisation == "remote" and download_obj.container_system == "singularity"
        assert isinstance(download_obj.containers_remote, list) and len(download_obj.containers_remote) == 0
        # read in the file
        containers_remote = SingularityFetcher.read_remote_singularity_containers(download_obj.container_cache_index)
        assert len(containers_remote) == 33
        assert "depot.galaxyproject.org-singularity-salmon-1.5.2--h84f40af_0.img" in containers_remote
        assert "MV Rena" not in containers_remote  # decoy in test file


#
# Tests for the FileDownloader class
#
class FileDownloaderTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    #
    # Test for 'download_file'
    #
    @with_temporary_folder
    def test_file_download(self, outdir):
        outdir = Path(outdir)
        with ContainerProgress() as progress:
            downloader = FileDownloader(progress)

            # Activate the caplog: all download attempts must be logged (even failed ones)
            self._caplog.clear()
            with self._caplog.at_level(logging.DEBUG):
                # No task initially
                assert progress.tasks == []
                assert progress._task_index == 0

                # Download a file
                src_url = "https://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.fasta.fai"
                output_path = outdir / Path(src_url).name
                downloader.download_file(src_url, output_path)
                assert (output_path).exists()
                assert os.path.getsize(output_path) == 27
                assert (
                    "nf_core.pipelines.download.singularity",
                    logging.DEBUG,
                    f"Downloading '{src_url}' to '{output_path}'",
                ) in self._caplog.record_tuples

                # A task was added but is now gone
                assert progress._task_index == 1
                assert progress.tasks == []

                # No content at the URL
                src_url = "http://www.google.com/generate_204"
                output_path = outdir / Path(src_url).name
                with pytest.raises(DownloadError):
                    downloader.download_file(src_url, output_path)
                assert not (output_path).exists()
                assert (
                    "nf_core.pipelines.download.singularity",
                    logging.DEBUG,
                    f"Downloading '{src_url}' to '{output_path}'",
                ) in self._caplog.record_tuples

                # A task was added but is now gone
                assert progress._task_index == 2
                assert progress.tasks == []

                # Invalid URL (schema)
                src_url = "dummy://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.fasta.fax"
                output_path = outdir / Path(src_url).name
                with pytest.raises(requests.exceptions.InvalidSchema):
                    downloader.download_file(src_url, output_path)
                assert not (output_path).exists()
                assert (
                    "nf_core.pipelines.download.singularity",
                    logging.DEBUG,
                    f"Downloading '{src_url}' to '{output_path}'",
                ) in self._caplog.record_tuples

                # A task was added but is now gone
                assert progress._task_index == 3
                assert progress.tasks == []

            # Fire in the hole ! The download will be aborted and no output file will be created
            src_url = "https://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.fasta.fai"
            output_path = outdir / Path(src_url).name
            os.unlink(output_path)
            downloader.kill_with_fire = True
            with pytest.raises(KeyboardInterrupt):
                downloader.download_file(src_url, output_path)
            assert not (output_path).exists()

    #
    # Test for 'download_files_in_parallel'
    #
    @with_temporary_folder
    def test_parallel_downloads(self, outdir):
        outdir = Path(outdir)

        # Prepare the download paths
        def make_tuple(url):
            return (url, (outdir / Path(url).name))

        download_fai = make_tuple(
            "https://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.fasta.fai"
        )
        download_dict = make_tuple(
            "https://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.dict"
        )
        download_204 = make_tuple("http://www.google.com/generate_204")
        download_schema = make_tuple(
            "dummy://github.com/nf-core/test-datasets/raw/refs/heads/modules/data/genomics/sarscov2/genome/genome.fasta.fax"
        )

        with ContainerProgress() as progress:
            downloader = FileDownloader(progress)

            # Download two files
            assert downloader.kill_with_fire is False
            downloads = [download_fai, download_dict]
            downloaded_files = downloader.download_files_in_parallel(downloads, parallel_downloads=1)
            assert len(downloaded_files) == 2
            assert downloaded_files == downloads
            assert (download_fai[1]).exists()
            assert (download_dict[1]).exists()
            assert downloader.kill_with_fire is False
            (download_fai[1]).unlink()
            (download_dict[1]).unlink()

            # This time, the second file will raise an exception
            assert downloader.kill_with_fire is False
            downloads = [download_fai, download_204]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1)
            assert downloader.kill_with_fire is False
            assert (download_fai[1]).exists()
            assert not (download_204[1]).exists()
            (download_fai[1]).unlink()

            # Now we swap the two files. The first one will raise an exception but the
            # second one will still be downloaded because only KeyboardInterrupt can
            # stop everything altogether.
            assert downloader.kill_with_fire is False
            downloads = [download_204, download_fai]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1)
            assert downloader.kill_with_fire is False
            assert (download_fai[1]).exists()
            assert not (download_204[1]).exists()
            (download_fai[1]).unlink()

            # We check that there's the same behaviour with `requests` errors.
            assert downloader.kill_with_fire is False
            downloads = [download_schema, download_fai]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1)
            assert downloader.kill_with_fire is False
            assert (download_fai[1]).exists()
            assert not (download_schema[1]).exists()
            (download_fai[1]).unlink()

            # Now we check the callback method
            callbacks = []

            def callback(*args):
                callbacks.append(args)

            # We check the same scenarios as above
            callbacks = []
            downloads = [download_fai, download_dict]
            downloader.download_files_in_parallel(downloads, parallel_downloads=1, callback=callback)
            assert len(callbacks) == 2
            assert callbacks == [
                (download_fai, FileDownloader.Status.DONE),
                (download_dict, FileDownloader.Status.DONE),
            ]

            callbacks = []
            downloads = [download_fai, download_204]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1, callback=callback)
            assert len(callbacks) == 2
            assert callbacks == [
                (download_fai, FileDownloader.Status.DONE),
                (download_204, FileDownloader.Status.ERROR),
            ]

            callbacks = []
            downloads = [download_204, download_fai]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1, callback=callback)
            assert len(callbacks) == 2
            assert callbacks == [
                (download_204, FileDownloader.Status.ERROR),
                (download_fai, FileDownloader.Status.DONE),
            ]

            callbacks = []
            downloads = [download_schema, download_fai]
            with pytest.raises(DownloadError):
                downloader.download_files_in_parallel(downloads, parallel_downloads=1, callback=callback)
            assert len(callbacks) == 2
            assert callbacks == [
                (download_schema, FileDownloader.Status.ERROR),
                (download_fai, FileDownloader.Status.DONE),
            ]

            # Finally, we check how the function behaves when a KeyboardInterrupt is raised
            with mock.patch("concurrent.futures.wait", side_effect=KeyboardInterrupt):
                callbacks = []
                downloads = [download_fai, download_204, download_dict]
                with pytest.raises(KeyboardInterrupt):
                    downloader.download_files_in_parallel(downloads, parallel_downloads=1, callback=callback)
                assert len(callbacks) == 3
                # Note: whn the KeyboardInterrupt is raised, download_204 and download_dict are not yet started.
                # They are therefore cancelled and pushed to the callback list immediately. download_fai is last
                # because it is running and can't be cancelled.
                assert callbacks == [
                    (download_204, FileDownloader.Status.CANCELLED),
                    (download_dict, FileDownloader.Status.CANCELLED),
                    (download_fai, FileDownloader.Status.ERROR),
                ]
