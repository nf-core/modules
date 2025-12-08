"""Downloads a nf-core pipeline to the local file system."""

import io
import json
import logging
import os
import re
import shutil
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Any, Literal
from zipfile import ZipFile

import questionary
import requests
import rich

import nf_core
import nf_core.pipelines.list
import nf_core.utils
from nf_core.pipelines.download.container_fetcher import ContainerFetcher
from nf_core.pipelines.download.docker import DockerFetcher
from nf_core.pipelines.download.singularity import SINGULARITY_CACHE_DIR_ENV_VAR, SingularityFetcher
from nf_core.pipelines.download.utils import DownloadError, intermediate_dir_with_cd
from nf_core.pipelines.download.workflow_repo import WorkflowRepo
from nf_core.utils import (
    NF_INSPECT_MIN_NF_VERSION,
    NFCORE_VER_LAST_WITHOUT_NF_INSPECT,
    check_nextflow_version,
    gh_api,
    pretty_nf_version,
    run_cmd,
)

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=nf_core.utils.rich_force_colors(),
)


class DownloadWorkflow:
    """Downloads a nf-core workflow from GitHub to the local file system.

    Can also download its Singularity container image if required.

    Args:
        pipeline (str | None): The name of an nf-core-compatible pipeline in the form org/repo
        revision (tuple[str] | str | None): The workflow revision(s) to download, like `1.0` or `dev` . Defaults to None.
        outdir (Path | None): Path to the local download directory. Defaults to None.
        compress_type (str | None): Type of compression for the downloaded files. Defaults to None.
        force (bool): Flag to force download even if files already exist (overwrite existing files). Defaults to False.
        platform (bool): Flag to customize the download for Seqera Platform (convert to git bare repo). Defaults to False.
        download_configuration (str | None): Download the configuration files from nf-core/configs. Defaults to None.
        additional_tags (list[str] | str | None): Specify additional tags to add to the downloaded pipeline. Defaults to None.
        container_system (str): The container system to use (e.g., "singularity"). Defaults to None.
        container_library (tuple[str] | str | None): The container libraries (registries) to use. Defaults to None.
        container_cache_utilisation (str | None): If a local or remote cache of already existing container images should be considered. Defaults to None.
        container_cache_index (Path | None): An index for the remote container cache. Defaults to None.
        parallel (int): The number of parallel downloads to use. Defaults to 4.
        hide_progress (bool): Flag to hide the progress bar. Defaults to False.
    """

    def __init__(
        self,
        pipeline: str | None = None,
        revision: tuple[str, ...] | str | None = None,
        outdir: Path | None = None,
        compress_type: str | None = None,
        force: bool = False,
        platform: bool = False,
        download_configuration: str | None = None,
        additional_tags: tuple[str, ...] | str | None = None,
        container_system: str | None = None,
        container_library: tuple[str, ...] | str | None = None,
        container_cache_utilisation: str | None = None,
        container_cache_index: Path | None = None,
        parallel: int = 4,
        hide_progress: bool = False,
    ):
        # Verify that the flags provided make sense together
        if (
            container_system == "docker"
            and container_cache_utilisation != "copy"
            and container_cache_utilisation is not None
        ):
            raise DownloadError(
                "Only the 'copy' option for --container-cache-utilisation is supported for Docker images. "
            )

        self._pipeline = pipeline
        if isinstance(revision, str):
            self.revision = [revision]
        elif isinstance(revision, tuple):
            self.revision = [*revision]
        else:
            self.revision = []
        self._outdir: Path | None = Path(outdir) if outdir is not None else None
        self.output_filename: Path | None = None

        self.compress_type = compress_type
        self.force = force
        self.hide_progress = hide_progress
        self.platform = platform
        self.fullname: str | None = None
        # downloading configs is not supported for Seqera Platform downloads.
        self.include_configs = True if download_configuration == "yes" and not bool(platform) else False
        # Additional tags to add to the downloaded pipeline. This enables to mark particular commits or revisions with
        # additional tags, e.g. "stable", "testing", "validated", "production" etc. Since this requires a git-repo, it is only
        # available for the bare / Seqera Platform download.
        self.additional_tags: list[str] | None
        if isinstance(additional_tags, str) and bool(len(additional_tags)) and self.platform:
            self.additional_tags = [additional_tags]
        elif isinstance(additional_tags, tuple) and bool(len(additional_tags)) and self.platform:
            self.additional_tags = [*additional_tags]
        else:
            self.additional_tags = None

        self.container_system = container_system
        self.container_fetcher: ContainerFetcher | None = None
        # Check if a cache or libraries were specfied even though singularity was not
        if self.container_system != "singularity":
            if container_cache_index:
                log.warning("The flag '--container-cache-index' is set, but not selected to fetch singularity images")
                self.prompt_use_singularity(
                    "The '--container-cache-index' flag is only applicable when fetching singularity images"
                )

            if container_library:
                log.warning("You have specified container libraries but not selected to fetch singularity image")
                self.prompt_use_singularity(
                    "The '--container-library' flag is only applicable when fetching singularity images"
                )  # Is this correct?

        # Manually specified container library (registry)
        if isinstance(container_library, str) and bool(len(container_library)):
            self.container_library = [container_library]
        elif isinstance(container_library, tuple) and bool(len(container_library)):
            self.container_library = [*container_library]
        else:
            self.container_library = ["quay.io"]
        # Create a new set and add all values from self.container_library (CLI arguments to --container-library)
        self.registry_set = set(self.container_library) if hasattr(self, "container_library") else set()
        # if a container_cache_index is given, use the file and overrule choice.
        self.container_cache_utilisation = "remote" if container_cache_index else container_cache_utilisation
        self.container_cache_index = container_cache_index
        # allows to specify a container library / registry or a respective mirror to download images from
        self.parallel = parallel
        self.hide_progress = hide_progress

        if not gh_api.has_init:
            gh_api.lazy_init()
        self.authenticated = gh_api.auth is not None

        self.wf_revisions: list[dict[str, Any]] = []
        self.wf_branches: dict[str, Any] = {}
        self.wf_sha: dict[str, str] = {}
        self.wf_download_url: dict[str, str] = {}
        self.nf_config: dict[str, str] = {}
        self.containers: list[str] = []
        self.containers_remote: list[str] = []  # stores the remote images provided in the file.

        # Fetch remote workflows
        self.wfs = nf_core.pipelines.list.Workflows()
        self.wfs.get_remote_workflows()

    @property
    def pipeline(self) -> str:
        """
        Get the pipeline name.
        """
        assert self._pipeline is not None  # mypy
        return self._pipeline

    @pipeline.setter
    def pipeline(self, pipeline: str) -> None:
        """
        Set the pipeline name.
        """
        self._pipeline = pipeline

    @property
    def outdir(self) -> Path:
        """
        Get the output directory for the download.
        """
        assert self._outdir is not None  # mypy
        return self._outdir

    @outdir.setter
    def outdir(self, outdir: Path) -> None:
        """
        Set the output directory for the download.
        """
        self._outdir = outdir

    def download_workflow(self) -> None:
        """Starts a nf-core workflow download."""

        # Get workflow details
        try:
            self.prompt_pipeline_name()
            self.pipeline, self.wf_revisions, self.wf_branches = nf_core.utils.get_repo_releases_branches(
                self.pipeline, self.wfs
            )
            self.prompt_revision()
            self.get_revision_hash()

            # After this point the outdir should be set
            assert self.outdir is not None  # mypy

            # Inclusion of configs is unnecessary for Seqera Platform.
            if not self.platform and self.include_configs is None:
                self.prompt_config_inclusion()
            # Prompt the user for whether containers should be downloaded
            if self.container_system is None:
                self.prompt_container_download()

            # Check if we have an outdated Nextflow version
            if (
                self.container_system is not None
                and self.container_system != "none"
                and not check_nextflow_version(NF_INSPECT_MIN_NF_VERSION)
            ):
                log.error(
                    f"Container download requires Nextflow version >= {pretty_nf_version(NF_INSPECT_MIN_NF_VERSION)}\n"
                    f"Please update your Nextflow version with [magenta]'nextflow self-update'[/]\n"
                    f"or use a version of 'nf-core/tools' <= {'.'.join([str(i) for i in NFCORE_VER_LAST_WITHOUT_NF_INSPECT])}"
                )
                return

            # Setup the appropriate ContainerFetcher object
            self.setup_container_fetcher()

            # Nothing meaningful to compress here.
            if not self.platform:
                self.prompt_compression_type()
        except AssertionError as e:
            raise DownloadError(e) from e

        summary_log = [
            f"Pipeline revision: '{', '.join(self.revision) if len(self.revision) < 5 else self.revision[0] + ',[' + str(len(self.revision) - 2) + ' more revisions],' + self.revision[-1]}'",
            f"Use containers: '{self.container_system}'",
        ]
        if self.container_system:
            summary_log.append(f"Container library: '{', '.join(self.container_library)}'")
        if self.container_system == "singularity" and os.environ.get(SINGULARITY_CACHE_DIR_ENV_VAR) is not None:
            summary_log.append(
                f"Using [blue]{SINGULARITY_CACHE_DIR_ENV_VAR}[/]': {os.environ[SINGULARITY_CACHE_DIR_ENV_VAR]}'"
            )
            if self.containers_remote:
                summary_log.append(
                    f"Successfully read {len(self.containers_remote)} containers from the remote '{SINGULARITY_CACHE_DIR_ENV_VAR}' contents."
                )

        # Set an output filename now that we have the outdir
        if self.platform:
            self.output_filename = self.outdir.parent / (self.outdir.name + ".git")
            summary_log.append(f"Output file: '{self.output_filename}'")
        elif self.compress_type is not None:
            self.output_filename = self.outdir.parent / (self.outdir.name + "." + self.compress_type)
            summary_log.append(f"Output file: '{self.output_filename}'")
        else:
            summary_log.append(f"Output directory: '{self.outdir}'")

        if not self.platform:
            # Only show entry, if option was prompted.
            summary_log.append(f"Include default institutional configuration: '{self.include_configs}'")
        else:
            summary_log.append(f"Enabled for Seqera Platform: '{self.platform}'")

        # Check that the outdir doesn't already exist
        if self.outdir.exists():
            if not self.force:
                raise DownloadError(
                    f"Output directory '{self.outdir}' already exists (use [red]--force[/] to overwrite)"
                )
            log.warning(f"Deleting existing output directory: '{self.outdir}'")
            shutil.rmtree(self.outdir)

        # Check that compressed output file doesn't already exist
        if self.output_filename and self.output_filename.exists():
            if not self.force:
                raise DownloadError(
                    f"Output file '{self.output_filename}' already exists (use [red]--force[/] to overwrite)"
                )
            log.warning(f"Deleting existing output file: '{self.output_filename}'")
            self.output_filename.unlink()

        # Summary log
        log_lines = "\n".join(summary_log)
        log.info(f"Saving '{self.pipeline}'\n{log_lines}")

        # Perform the actual download
        if self.platform:
            log.info("Downloading workflow for Seqera Platform")
            self.download_workflow_platform()
        else:
            log.info("Downloading workflow")
            self.download_workflow_static()

        # The container fetcher might have some clean-up code, call it
        if self.container_fetcher:
            self.container_fetcher.cleanup()

    def download_workflow_static(self) -> None:
        """Downloads a nf-core workflow from GitHub to the local file system in a self-contained manner."""

        # Download the centralised configs first
        if self.include_configs:
            log.info("Downloading centralised configs from GitHub")
            self.download_configs()

        # Download the pipeline files for each selected revision
        log.info("Downloading workflow files from GitHub")

        for revision, wf_sha, download_url in zip(self.revision, self.wf_sha.values(), self.wf_download_url.values()):
            revision_dirname = self.download_wf_files(revision=revision, wf_sha=wf_sha, download_url=download_url)

            if self.include_configs:
                try:
                    self.wf_use_local_configs(revision_dirname)
                except FileNotFoundError as e:
                    raise DownloadError("Error editing pipeline config file to use local configs!") from e

            # Collect all required container images
            if self.container_system in {"singularity", "docker"}:
                workflow_directory = self.outdir / revision_dirname
                self.find_container_images(workflow_directory, revision)

                try:
                    self.download_container_images(workflow_directory, revision)
                except OSError as e:
                    raise DownloadError(f"[red]{e}[/]") from e

        # Compress into an archive
        if self.compress_type is not None:
            log.info("Compressing output into archive")
            self.compress_download()

    def download_workflow_platform(self, location: Path | None = None) -> None:
        """Create a bare-cloned git repository of the workflow, so it can be launched with `tw launch` as file:/ pipeline"""
        assert self.output_filename is not None  # mypy

        log.info("Collecting workflow from GitHub")

        self.workflow_repo = WorkflowRepo(
            remote_url=f"https://github.com/{self.pipeline}.git",
            revision=self.revision if self.revision else None,
            commit=self.wf_sha.values() if bool(self.wf_sha) else None,
            additional_tags=self.additional_tags,
            location=location if location else None,  # manual location is required for the tests to work
            in_cache=False,
        )

        # Remove tags for those revisions that had not been selected
        self.workflow_repo.tidy_tags_and_branches()

        # create a bare clone of the modified repository needed for Seqera Platform
        self.workflow_repo.bare_clone(self.outdir / self.output_filename)

        # extract the required containers
        if self.container_system in {"singularity", "docker"}:
            for revision, commit in self.wf_sha.items():
                # Checkout the repo in the current revision
                self.workflow_repo.checkout(commit)
                # Collect all required singularity images
                workflow_directory = self.workflow_repo.access()
                self.find_container_images(workflow_directory, revision)

                try:
                    self.download_container_images(workflow_directory, revision)
                except OSError as e:
                    raise DownloadError(f"[red]{e}[/]") from e

        # Justify why compression is skipped for Seqera Platform downloads (Prompt is not shown, but CLI argument could have been set)
        if self.compress_type is not None:
            log.info(
                "Compression choice is ignored for Seqera Platform downloads since nothing can be reasonably compressed."
            )

    def prompt_pipeline_name(self) -> None:
        """Prompt for the pipeline name if not set with a flag"""

        if self._pipeline is None:
            stderr.print("Specify the name of a nf-core pipeline or a GitHub repository name (user/repo).")
            self.pipeline = nf_core.utils.prompt_remote_pipeline_name(self.wfs)

    def prompt_revision(self) -> None:
        """
        Prompt for pipeline revision / branch
        Prompt user for revision tag if '--revision' was not set
        If --platform is specified, allow to select multiple revisions
        Also the static download allows for multiple revisions, but
        we do not prompt this option interactively.
        """
        if not bool(self.revision):
            (choice, tag_set) = nf_core.utils.prompt_pipeline_release_branch(
                self.wf_revisions, self.wf_branches, multiple=self.platform
            )
            """
            The checkbox() prompt unfortunately does not support passing a Validator,
            so a user who keeps pressing Enter will flounder past the selection without choice.

            bool(choice), bool(tag_set):
            #############################
            True,  True:  A choice was made and revisions were available.
            False, True:  No selection was made, but revisions were available -> defaults to all available.
            False, False: No selection was made because no revisions were available -> raise AssertionError.
            True,  False: Congratulations, you found a bug! That combo shouldn't happen.
            """

            if bool(choice):
                # have to make sure that self.revision is a list of strings, regardless if choice is str or list of strings.
                (self.revision.append(choice) if isinstance(choice, str) else self.revision.extend(choice))
            else:
                if bool(tag_set):
                    self.revision = tag_set
                    log.info("No particular revision was selected, all available will be downloaded.")
                else:
                    raise AssertionError(f"No revisions of {self.pipeline} available for download.")

    def get_revision_hash(self) -> None:
        """Find specified revision / branch / commit hash"""

        for revision in self.revision:  # revision is a list of strings, but may be of length 1
            # Branch
            if revision in self.wf_branches.keys():
                self.wf_sha = {**self.wf_sha, revision: self.wf_branches[revision]}

            else:
                # Revision
                for r in self.wf_revisions:
                    if r["tag_name"] == revision:
                        self.wf_sha = {**self.wf_sha, revision: r["tag_sha"]}
                        break

                else:
                    # Commit - full or short hash
                    if commit_id := nf_core.utils.get_repo_commit(self.pipeline, revision):
                        self.wf_sha = {**self.wf_sha, revision: commit_id}
                        continue

                    # Can't find the revisions or branch - throw an error
                    log.info(
                        "Available {} revisions: '{}'".format(
                            self.pipeline,
                            "', '".join([r["tag_name"] for r in self.wf_revisions]),
                        )
                    )
                    log.info("Available {} branches: '{}'".format(self.pipeline, "', '".join(self.wf_branches.keys())))
                    raise AssertionError(
                        f"Not able to find revision / branch / commit '{revision}' for {self.pipeline}"
                    )

        # Set the outdir
        if not self._outdir:
            if len(self.wf_sha) > 1:
                self.outdir = Path(
                    f"{self.pipeline.replace('/', '-').lower()}_{datetime.now().strftime('%Y-%m-%d_%H-%M')}"
                )
            else:
                self.outdir = Path(f"{self.pipeline.replace('/', '-').lower()}_{self.revision[0]}")

        if not self.platform:
            for revision, wf_sha in self.wf_sha.items():
                # Set the download URL - only applicable for classic downloads
                if self.authenticated:
                    # For authenticated downloads, use the GitHub API
                    self.wf_download_url = {
                        **self.wf_download_url,
                        revision: f"https://api.github.com/repos/{self.pipeline}/zipball/{wf_sha}",
                    }
                else:
                    # For unauthenticated downloads, use the archive URL
                    self.wf_download_url = {
                        **self.wf_download_url,
                        revision: f"https://github.com/{self.pipeline}/archive/{wf_sha}.zip",
                    }

    def prompt_config_inclusion(self) -> None:
        """Prompt for inclusion of institutional configurations"""
        if stderr.is_interactive:  # Use rich auto-detection of interactive shells
            self.include_configs = questionary.confirm(
                "Include the nf-core's default institutional configuration files into the download?",
                style=nf_core.utils.nfcore_question_style,
            ).ask()
        else:
            self.include_configs = False
            # do not include by default.

    def prompt_container_download(self) -> None:
        """Prompt whether to download container images or not"""

        if self.container_system is None and stderr.is_interactive and not self.platform:
            stderr.print("\nIn addition to the pipeline code, this tool can download software containers.")
            self.container_system = questionary.select(
                "Download software container images:",
                choices=["none", "singularity", "docker"],
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

    def setup_container_fetcher(self) -> None:
        """
        Create the appropriate ContainerFetcher object
        """
        assert self.outdir is not None  # mypy
        try:
            if self.container_system == "singularity":
                self.container_fetcher = SingularityFetcher(
                    outdir=self.outdir,
                    container_library=self.container_library,
                    registry_set=self.registry_set,
                    container_cache_utilisation=self.container_cache_utilisation,
                    container_cache_index=self.container_cache_index,
                    parallel=self.parallel,
                    hide_progress=self.hide_progress,
                )
            elif self.container_system == "docker":
                self.container_fetcher = DockerFetcher(
                    outdir=self.outdir,
                    registry_set=self.registry_set,
                    container_library=self.container_library,
                    parallel=self.parallel,
                    hide_progress=self.hide_progress,
                )
            else:
                self.container_fetcher = None
        except OSError as e:
            raise DownloadError(e)

    def prompt_use_singularity(self, fail_message: str) -> None:
        use_singularity = questionary.confirm(
            "Do you want to download singularity images?",
            style=nf_core.utils.nfcore_question_style,
        ).ask()
        if use_singularity:
            self.container_system = "singularity"
        else:
            raise DownloadError(fail_message)

    def prompt_compression_type(self) -> None:
        """Ask user if we should compress the downloaded files"""
        if self.compress_type is None:
            stderr.print(
                "\nIf transferring the downloaded files to another system, it can be convenient to have everything compressed in a single file."
            )
            if self.container_system == "singularity":
                stderr.print(
                    "[bold]This is [italic]not[/] recommended when downloading Singularity images, as it can take a long time and saves very little space."
                )
            self.compress_type = questionary.select(
                "Choose compression type:",
                choices=[
                    "none",
                    "tar.gz",
                    "tar.bz2",
                    "zip",
                ],
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

        # Correct type for no-compression
        if self.compress_type == "none":
            self.compress_type = None

    def download_wf_files(self, revision: str, wf_sha: str, download_url: str) -> str:
        """Downloads workflow files from GitHub to the :attr:`self.outdir`."""

        log.debug(f"Downloading {download_url}")

        # GitHub API download: fetch via API and get topdir from zip contents
        content = gh_api.get(download_url).content
        with ZipFile(io.BytesIO(content)) as zipfile:
            topdir = zipfile.namelist()[0]  # API zipballs have a generated directory name
            zipfile.extractall(self.outdir)

        # Create a filesystem-safe version of the revision name for the directory
        revision_dirname = re.sub("[^0-9a-zA-Z]+", "_", revision)
        # Account for name collisions, if there is a branch / release named "configs" or container output dir
        if revision_dirname in ["configs", self.get_container_output_dir()]:
            revision_dirname = re.sub("[^0-9a-zA-Z]+", "_", self.pipeline + revision_dirname)

        # Rename the internal directory name to be more friendly
        (self.outdir / topdir).rename(self.outdir / revision_dirname)

        # Make downloaded files executable
        for dirpath, _, filelist in os.walk(self.outdir / revision_dirname):
            for fname in filelist:
                (Path(dirpath) / fname).chmod(0o775)

        return revision_dirname

    def download_configs(self) -> None:
        """Downloads the centralised config profiles from nf-core/configs to :attr:`self.outdir`."""
        configs_zip_url = "https://github.com/nf-core/configs/archive/master.zip"
        configs_local_dir = "configs-master"
        log.debug(f"Downloading {configs_zip_url}")

        # Download GitHub zip file into memory and extract
        url = requests.get(configs_zip_url)
        with ZipFile(io.BytesIO(url.content)) as zipfile:
            zipfile.extractall(self.outdir)

        # Rename the internal directory name to be more friendly
        (self.outdir / configs_local_dir).rename(self.outdir / "configs")

        # Make downloaded files executable

        for dirpath, _, filelist in os.walk(self.outdir / "configs"):
            for fname in filelist:
                (Path(dirpath) / fname).chmod(0o775)

    def wf_use_local_configs(self, revision_dirname: str) -> None:
        """Edit the downloaded nextflow.config file to use the local config files"""

        assert self.outdir is not None  # mypy
        nfconfig_fn = (self.outdir / revision_dirname) / "nextflow.config"
        find_str = "https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}"
        repl_str = "${projectDir}/../configs/"
        log.debug(f"Editing 'params.custom_config_base' in '{nfconfig_fn}'")

        # Load the nextflow.config file into memory
        with open(nfconfig_fn) as nfconfig_fh:
            nfconfig = nfconfig_fh.read()

        # Replace the target string
        log.debug(f"Replacing '{find_str}' with '{repl_str}'")
        nfconfig = nfconfig.replace(find_str, repl_str)

        # Append the singularity.cacheDir to the end if we need it
        if self.container_system == "singularity" and self.container_cache_utilisation == "copy":
            nfconfig += (
                f"\n\n// Added by `nf-core pipelines download` v{nf_core.__version__} //\n"
                + 'singularity.cacheDir = "${projectDir}/../singularity-images/"'
                + "\n///////////////////////////////////////"
            )

        # Write the file out again
        log.debug(f"Updating '{nfconfig_fn}'")
        with open(nfconfig_fn, "w") as nfconfig_fh:
            nfconfig_fh.write(nfconfig)

    def find_container_images(
        self, workflow_directory: Path, revision: str, with_test_containers: bool = True, entrypoint: str = "main.nf"
    ) -> None:
        """
        Find container image names for workflow using the `nextflow inspect` command.

        Requires Nextflow >= 25.04.4

        Args:
            workflow_directory (Path): The directory containing the workflow files.
            entrypoint (str): The entrypoint for the `nextflow inspect` command.
        """

        log.info(
            f"Fetching container names for workflow revision {revision} using [magenta bold]nextflow inspect[/]. This might take a while."
        )
        try:
            # TODO: Select container system via profile. Is this stable enough?
            # NOTE: We will likely don't need this after the switch to Seqera containers
            profile_str = f"{self.container_system}"
            if with_test_containers:
                profile_str += ",test,test_full"
            profile = f"-profile {profile_str}" if self.container_system else ""

            working_dir = Path().absolute()
            with intermediate_dir_with_cd(working_dir):
                # Run nextflow inspect
                executable = "nextflow"
                cmd_params = f"inspect -format json {profile} {working_dir / workflow_directory / entrypoint}"
                cmd_out = run_cmd(executable, cmd_params)
                if cmd_out is None:
                    raise DownloadError("Failed to run `nextflow inspect`. Please check your Nextflow installation.")

                out, _ = cmd_out
                out_json = json.loads(out)
                # NOTE: Should we save the container name too to have more meta information?
                named_containers = {proc["name"]: proc["container"] for proc in out_json["processes"]}
                # We only want to process unique containers
                self.containers = list(set(named_containers.values()))

        except RuntimeError as e:
            log.error("Running 'nextflow inspect' failed with the following error")
            raise DownloadError(e)

        except KeyError as e:
            log.error("Failed to parse output of 'nextflow inspect' to extract containers")
            raise DownloadError(e)

    def get_container_output_dir(self) -> Path:
        assert self.outdir is not None  # mypy
        return self.outdir / f"{self.container_system}-images"

    def download_container_images(self, workflow_directory: Path, current_revision: str = "") -> None:
        """
        Fetch the container images with the appropriate ContainerFetcher

        Args:
            current_revision (str): The current revision of the workflow.
        """

        if len(self.containers) == 0:
            log.info("No container names found in workflow")
        else:
            log.info(
                f"Processing workflow revision {current_revision}, found {len(self.containers)} container image{'s' if len(self.containers) > 1 else ''} in total."
            )
            log.debug(f"Container names: {self.containers}")

            out_path_dir = self.get_container_output_dir().absolute()

            # Check that the directories exist
            if not out_path_dir.is_dir():
                log.debug(f"Output directory not found, creating: {out_path_dir}")
                out_path_dir.mkdir(parents=True)

            if self.container_fetcher is not None:
                self.container_fetcher.fetch_containers(self.containers, self.containers_remote, workflow_directory)

    def compress_download(self) -> None:
        """Take the downloaded files and make a compressed .tar.gz archive."""
        log.debug(f"Creating archive: {self.output_filename}")

        # .tar.gz and .tar.bz2 files
        if self.compress_type in ["tar.gz", "tar.bz2"]:
            ctype_and_mode: dict[str, tuple[Literal["gz", "bz2"], Literal["w:gz", "w:bz2"]]] = {
                "tar.gz": ("gz", "w:gz"),
                "tar.bz2": ("bz2", "w:bz2"),
            }  # This ugly thing is required for typing
            ctype, mode = ctype_and_mode[self.compress_type]
            with tarfile.open(self.output_filename, mode) as tar:
                tar.add(self.outdir, arcname=self.outdir.name)
            tar_flags = "xzf" if ctype == "gz" else "xjf"
            log.info(f"Command to extract files: [bright_magenta]tar -{tar_flags} {self.output_filename}[/]")

        # .zip files
        if self.compress_type == "zip":
            assert self.output_filename is not None  # mypy
            with ZipFile(self.output_filename, "w") as zip_file:
                # Iterate over all the files in directory
                for folder_name, _, filenames in os.walk(self.outdir):
                    for filename in filenames:
                        # create complete filepath of file in directory
                        file_path = Path(folder_name) / filename
                        # Add file to zip
                        zip_file.write(file_path)
            log.info(f"Command to extract files: [bright_magenta]unzip {self.output_filename}[/]")

        # Delete original files
        log.debug(f"Deleting uncompressed files: '{self.outdir}'")
        shutil.rmtree(self.outdir)

        # Calculate md5sum for output file
        log.info(f"MD5 checksum for '{self.output_filename}': [blue]{nf_core.utils.file_md5(self.output_filename)}[/]")
