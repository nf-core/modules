import contextlib
import logging
import re
import shutil
from abc import ABC, abstractmethod
from collections.abc import Callable, Collection, Container, Generator, Iterable
from pathlib import Path

import rich.progress

import nf_core.utils
from nf_core.pipelines.download.utils import intermediate_file

log = logging.getLogger(__name__)


class ContainerProgress(rich.progress.Progress):
    """
    Custom Progress bar class, allowing us to have two progress
    bars with different columns / layouts.
    Also provide helper functions to control the top-level task.
    """

    main_task: rich.progress.TaskID | None = None
    remote_fetch_task: rich.progress.TaskID | None = None
    remote_fetch_task_containers: list[str] | None = []
    copy_task: rich.progress.TaskID | None = None
    copy_task_containers: list[str] | None = []

    def __init__(self, disable=False):
        super().__init__(disable=disable)

    def get_task_types_and_columns(self):
        """
        Gets the possible task types for the progress bar.
        """
        task_types_and_columns = {
            "summary": (
                "[magenta]{task.description}",
                rich.progress.BarColumn(bar_width=None),
                "[progress.percentage]{task.percentage:>3.0f}%",
                "•",
                "[green]{task.completed}/{task.total} tasks completed",
            ),
            "remote_fetch": (
                "[cyan]{task.description}",
                rich.progress.BarColumn(bar_width=None),
                "[progress.percentage]{task.percentage:>3.0f}%",
                "•",
                "[green]{task.completed}/{task.total} tasks completed",
            ),
            "copy": (
                "[steel_blue]{task.description}",
                rich.progress.BarColumn(bar_width=None),
                "[progress.percentage]{task.percentage:>3.0f}%",
                "•",
                "[green]{task.completed}/{task.total} tasks completed",
            ),
        }
        return task_types_and_columns

    def get_renderables(self) -> Generator[rich.table.Table, None, None]:
        self.columns: Iterable[str | rich.progress.ProgressColumn]
        for task in self.tasks:
            for task_type, columns in self.get_task_types_and_columns().items():
                if task.fields.get("progress_type") == task_type:
                    self.columns = columns

            yield self.make_tasks_table([task])

    # These two functions allow callers not having to track the main TaskID
    # They are pass-through functions to the rich.progress methods
    def add_main_task(self, **kwargs) -> rich.progress.TaskID:
        """
        Add a top-level task to the progress bar.
        This task will be used to track the overall progress of the container downloads.
        """
        self.main_task = self.add_task(
            progress_type="summary",
            description="Processing container images",
            **kwargs,
        )
        return self.main_task

    def update_main_task(self, **kwargs) -> None:
        """
        Update the top-level task with new information.
        """
        self.update(self.main_task, **kwargs)

    def remove_main_task(self) -> None:
        """
        Remove the top-level task
        """
        self.remove_task(self.main_task)

    def add_remote_fetch_task(self, total: int, **kwargs) -> rich.progress.TaskID:
        """
        Add a task to the progress bar to track the progress of fetching remote containers.
        """
        self.remote_fetch_task = self.add_task(
            progress_type="remote_fetch",
            total=total,
            completed=0,
            **kwargs,
        )
        return self.remote_fetch_task

    def advance_remote_fetch_task(self) -> None:
        """
        Advance the remote fetch task, and if the container should not
        be copied then also advance the main task.
        """
        self.update(self.remote_fetch_task, advance=1)
        self.update_main_task(advance=1)

    def update_remote_fetch_task(self, **kwargs) -> None:
        """
        Update the remote fetch task with new information
        """
        self.update(self.remote_fetch_task, **kwargs)

    def remove_remote_fetch_task(self) -> None:
        """
        Remove the remote fetch task
        """
        self.remove_task(self.remote_fetch_task)

    def add_copy_task(self, total: int, **kwargs) -> rich.progress.TaskID:
        """
        Add a task to the progress bar to track the progress of copying containers.
        """
        self.copy_task = self.add_task(
            progress_type="copy",
            total=total,
            completed=0,
            **kwargs,
        )
        return self.copy_task

    def advance_copy_task(self) -> None:
        """
        Advance the copy task, and with it the main task.
        """
        self.update(self.copy_task, advance=1)
        self.update_main_task(advance=1)

    def update_copy_task(self, **kwargs) -> None:
        """
        Update the remote fetch task with new information
        """
        self.update(self.copy_task, **kwargs)

    def remove_copy_task(self) -> None:
        """
        Remove the copy task
        """
        self.remove_task(self.copy_task)

    @contextlib.contextmanager
    def sub_task(self, *args, **kwargs) -> Generator[rich.progress.TaskID, None, None]:
        """
        Context manager to create a sub-task under the main task.
        """
        task = self.add_task(*args, **kwargs)
        try:
            yield task
        finally:
            self.remove_task(task)


class ContainerFetcher(ABC):
    """
    Abstract class to manage all operations for fetching containers.

    It is currently subclasses by the SingularityFetcher and DockerFetcher classes,
    for fetching Singularity and Docker containers respectively.

    The guiding principles are that:
      - Container download/pull/copy methods are unaware of the concepts of
        "library" and "cache". They are just told to fetch a container and
        put it in a certain location.
      - Only the `fetch_containers` method is aware of the concepts of "library"
        and "cache". It is a sort of orchestrator that decides where to fetch
        each container and calls the appropriate methods.
      - All methods are integrated with a progress bar

    Args:
        container_output_dir (Path): The final destination for the container images.
        container_library (Iterable[str]): A collection of container libraries to use
        registry_set (Iterable[str]): A collection of registries to consider
        progress_factory (Callable[[], ContainerProgress]): A factory to create a progress bar.
        library_dir (Path | None): The directory to look for container images in.
        cache_dir (Path | None): A directory where container images might be cached.
        amend_cachedir (bool): Whether to amend the cache directory with the container images.
        parallel (int): The number of containers to fetch in parallel.
    """

    def __init__(
        self,
        container_output_dir: Path,
        container_library: Iterable[str],
        registry_set: Iterable[str],
        progress_factory: Callable[[bool], ContainerProgress],
        library_dir: Path | None,
        cache_dir: Path | None,
        amend_cachedir: bool,
        parallel: int = 4,
        hide_progress: bool = False,
    ) -> None:
        self._container_output_dir = container_output_dir
        self.container_library = list(container_library)
        self.base_registry_set: set[str] = set(registry_set)
        self._registry_set: set[str] | None = None

        self.kill_with_fire = False
        self.implementation: str | None = None
        self.name = None
        self.library_dir = library_dir
        self.cache_dir = cache_dir
        self.amend_cachedir = amend_cachedir
        self.parallel = parallel

        self.hide_progress = hide_progress
        self.progress_factory = progress_factory
        self.progress: ContainerProgress | None = None

    @property
    def progress(self) -> rich.progress.Progress:
        assert self._progress is not None  # mypy
        return self._progress

    @progress.setter
    def progress(self, progress: ContainerProgress | None) -> None:
        self._progress = progress

    @property
    def registry_set(self) -> set[str]:
        """
        Get the set of registries to use for the container download
        """
        assert self._registry_set is not None  # mypy
        return self._registry_set

    @registry_set.setter
    def registry_set(self, registry_set: set[str]) -> None:
        self._registry_set = registry_set

    def get_container_output_dir(self) -> Path:
        """
        Get the output directory for the container images.
        """
        return self._container_output_dir

    @abstractmethod
    def check_and_set_implementation(self) -> None:
        """
        Check if the container system is installed and available.

        Should update the `self.implementation` attribute with the found implementation

        Raises:
            OSError: If the container system is not installed or not in $PATH.
        """
        pass

    @abstractmethod
    def gather_registries(self, workflow_directory: Path) -> set[str]:
        """
        Gather the registries from the pipeline config and CLI arguments and store them in a set.

        Returns:
            set[str]: The set of registries.
        """
        pass

    def gather_config_registries(self, workflow_directory: Path, registry_keys: list[str] = []) -> set[str]:
        """
        Gather the registries from the pipeline config and store them in a set.

        Args:
            workflow_directory (Path): The directory containing the pipeline files we are currently processing
            registry_keys (list[str]): The list of registry keys to fetch from the pipeline config

        Returns:
            set[str]: The set of registries defined in the pipeline config
        """
        # Fetch the pipeline config
        nf_config = nf_core.utils.fetch_wf_config(workflow_directory)

        config_registries = set()
        for registry_key in registry_keys:
            if registry_key in nf_config:
                config_registries.add(nf_config[registry_key])

        return config_registries

    @abstractmethod
    def clean_container_file_extension(self, container_fn: str) -> str:
        """
        Clean the file extension of a container filename.

        Example implementation:

            # Detect file extension
            extension = ".img"
            if ".sif:" in out_name:
                extension = ".sif"
                out_name = out_name.replace(".sif:", "-")
            elif out_name.endswith(".sif"):
                extension = ".sif"
                out_name = out_name[:-4]
            # Strip : and / characters
            out_name = out_name.replace("/", "-").replace(":", "-")
            # Add file extension
            out_name = out_name + extension

        Args:
            container_fn (str): The filename of the container.
        Returns:
            str: The cleaned filename with the appropriate extension.
        """
        pass

    # We have dropped the explicit registries from the modules in favor of the configurable registries.
    # Unfortunately, Nextflow still expects the registry to be part of the file name, so we need functions
    # to support accessing container images with different registries (or no registry).
    def get_container_filename(self, container: str) -> str:
        """Return the expected filename for a container.

        Supports docker, http, oras, and singularity URIs in `container`.

        Registry names provided in `registries` are removed from the filename to ensure that the same image
        is used regardless of the registry. Only registry names that are part of `registries` are considered.
        If the image name contains another registry, it will be kept in the filename.

        For instance, docker.io/nf-core/ubuntu:20.04 will be nf-core-ubuntu-20.04.img *only* if the registry
        contains "docker.io".
        """

        # Generate file paths
        # Based on simpleName() function in Nextflow code:
        # https://github.com/nextflow-io/nextflow/blob/671ae6d85df44f906747c16f6d73208dbc402d49/modules/nextflow/src/main/groovy/nextflow/container/SingularityCache.groovy#L69-L94
        out_name = container
        # Strip URI prefix
        out_name = re.sub(r"^.*:\/\/", "", out_name)

        # Clean the file extension. This method must be implemented
        # by any subclass
        out_name = self.clean_container_file_extension(out_name)

        # Trim potential registries from the name for consistency.
        # This will allow pipelines to work offline without symlinked images,
        # if docker.registry / singularity.registry are set to empty strings at runtime, which can be included in the HPC config profiles easily.
        if self.registry_set:
            # Create a regex pattern from the set of registries
            trim_pattern = "|".join(f"^{re.escape(registry)}-?".replace("/", "[/-]") for registry in self.registry_set)
            # Use the pattern to trim the string
            out_name = re.sub(f"{trim_pattern}", "", out_name)

        return out_name

    def fetch_containers(
        self,
        containers: Collection[str],
        exclude_list: Container[str],
        workflow_directory: Path,
    ):
        """
        This is the main entrypoint of the container fetcher. It goes through
        all the containers we find and does the appropriate action; copying
        from cache or fetching from a remote location
        """

        # Create a new progress bar
        self.progress = self.progress_factory(self.hide_progress)

        # Collect registries defined in the workflow directory
        self.registry_set = self.gather_registries(workflow_directory)

        with self.progress:
            # Check each container in the list and defer actions
            containers_remote_fetch: list[tuple[str, Path]] = []
            containers_copy: list[tuple[str, Path, Path]] = []

            # The first task is to check what to do with each container
            total_tasks = len(containers)
            self.progress.add_main_task(total=total_tasks)

            for container in containers:
                container_filename = self.get_container_filename(container)

                # Files in the remote cache are already downloaded and can be ignored
                if container_filename in exclude_list:
                    log.debug(f"Skipping download of container '{container_filename}' as it is cached remotely.")
                    self.progress.update_main_task(advance=1)
                    continue

                # Generate file paths for all three locations
                output_path = self.get_container_output_dir() / container_filename

                if output_path.exists():
                    log.debug(
                        f"Skipping download of container '{container_filename}' as it is already in `{self.get_container_output_dir()}`."
                    )
                    self.progress.update_main_task(advance=1)
                    continue

                library_path = self.library_dir / container_filename if self.library_dir is not None else None
                cache_path = self.cache_dir / container_filename if self.cache_dir is not None else None

                # Get the container from the library
                if library_path and library_path.exists():
                    # Update the cache if needed
                    if cache_path and not cache_path.exists() and self.amend_cachedir:
                        containers_copy.append((container, library_path, cache_path))

                    if not self.amend_cachedir:
                        # We are not just amending the cache directory, so the file should be copied to the output
                        containers_copy.append((container, library_path, output_path))

                # Get the container from the cache
                elif cache_path and cache_path.exists() and not self.amend_cachedir:
                    log.debug(f"Container '{container_filename}' found in cache at '{cache_path}'.")
                    containers_copy.append((container, cache_path, output_path))

                # Image is not in library or cache
                else:
                    # Fetching of remote containers, either pulling or downloading, differs between docker and singularity:
                    # - Singularity images can either be downloaded from an http address, or pulled from a registry with `(singularity|apptainer) pull`
                    # - Docker images are always pulled, but needs the additional `docker image save` command for the image to be saved in the correct place
                    if cache_path:
                        # Download into the cache
                        containers_remote_fetch.append((container, cache_path))

                        # Do not copy to the output directory "(docker|singularity)-images" if we are solely amending the cache
                        if not self.amend_cachedir:
                            containers_copy.append((container, cache_path, output_path))
                            total_tasks += 1
                    else:
                        # There is no cache directory so download or pull directly to the output
                        containers_remote_fetch.append((container, output_path))

            self.progress.update_main_task(total=total_tasks)

            # Fetch containers from a remote location
            if containers_remote_fetch:
                self.progress.add_remote_fetch_task(
                    total=len(containers_remote_fetch),
                    description=f"Fetch remote {self.implementation} images",
                )
                self.fetch_remote_containers(containers_remote_fetch, parallel=self.parallel)
                self.progress.remove_remote_fetch_task()

            # Copy containers
            if containers_copy:
                self.progress.add_copy_task(
                    total=len(containers_copy),
                    description="Copy container images from/to cache",
                )
                for container, src_path, dest_path in containers_copy:
                    self.copy_image(container, src_path, dest_path)
                    self.progress.advance_copy_task()
                self.progress.remove_copy_task()

            self.progress.remove_main_task()
        # Unset the progress bar, so that we get an AssertionError if we access it after it is closed
        self.progress = None

    @abstractmethod
    def fetch_remote_containers(self, containers: list[tuple[str, Path]], parallel: int = 4) -> None:
        """
        Fetch remote containers

        - Singularity: pull or download images, depending on what address we have
        - Docker: pull and save images

        This function should update the main progress task accordingly
        """
        pass

    def copy_image(self, container: str, src_path: Path, dest_path: Path) -> None:
        """Copy container image from one directory to another."""
        # Check that the source path exists
        if not src_path.exists():
            log.error(f"Image '{container}' does not exist")
            return

        with intermediate_file(dest_path) as dest_path_tmp:
            shutil.copyfile(src_path, dest_path_tmp.name)

    def cleanup(self) -> None:
        """
        Cleanup any temporary files or resources.
        """
        pass
