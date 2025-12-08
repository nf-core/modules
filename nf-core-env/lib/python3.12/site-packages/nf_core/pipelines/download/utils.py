import contextlib
import importlib.resources
import logging
import os
import shutil
import tempfile
from collections.abc import Generator
from enum import Enum
from pathlib import Path

log = logging.getLogger(__name__)


class ContainerRegistryUrls(Enum):
    SEQERA_DOCKER = "community.wave.seqera.io/library"
    SEQERA_SINGULARITY = "community-cr-prod.seqera.io/docker/registry/v2"
    GALAXY_SINGULARITY = "depot.galaxyproject.org/singularity"


def copy_container_load_scripts(container_system: str, dest_dir: Path, make_exec: bool = True) -> tuple[str, Path]:
    container_load_scripts_subpackage = "nf_core.pipelines.download.load_scripts"
    script_name = f"{container_system}-load.sh"
    dest_path = dest_dir / script_name
    with importlib.resources.open_text(container_load_scripts_subpackage, script_name) as src:
        with open(dest_path, "w") as dest:
            shutil.copyfileobj(src, dest)
    if make_exec:
        dest_path.chmod(0o775)
    return script_name, dest_path


class DownloadError(RuntimeError):
    """A custom exception that is raised when nf-core pipelines download encounters a problem that we already took into consideration.
    In this case, we do not want to print the traceback, but give the user some concise, helpful feedback instead.
    """


@contextlib.contextmanager
def intermediate_file(output_path: Path) -> Generator[tempfile._TemporaryFileWrapper, None, None]:
    """Context manager to help ensure the output file is either complete or non-existent.
    It does that by creating a temporary file in the same directory as the output file,
    letting the caller write to it, and then moving it to the final location.
    If an exception is raised, the temporary file is deleted and the output file is not touched.
    """
    if output_path.is_dir():
        raise DownloadError(f"Output path '{output_path}' is a directory")
    if output_path.is_symlink():
        raise DownloadError(f"Output path '{output_path}' is a symbolic link")

    tmp = tempfile.NamedTemporaryFile(dir=output_path.parent, delete=False)
    try:
        yield tmp
        tmp.close()
        Path(tmp.name).rename(output_path)
    except:
        tmp_path = Path(tmp.name)
        if tmp_path.exists():
            tmp_path.unlink()
        raise


@contextlib.contextmanager
def intermediate_file_no_creation(output_path: Path) -> Generator[Path, None, None]:
    """
    Context manager to help ensure the output file is either complete or non-existent.

    'singularity/apptainer pull' requires that the output file does not exist before it is run.
    For pulling container we therefore create a temporary directory with and write to a file named
    'tempfile' in it. If the pull command is successful, we rename the temporary file to the output path.
    """
    if output_path.is_dir():
        raise DownloadError(f"Output path '{output_path}' is a directory")
    if output_path.is_symlink():
        raise DownloadError(f"Output path '{output_path}' is a symbolic link")

    tmp = tempfile.TemporaryDirectory(dir=output_path.parent)
    tmp_fn = Path(tmp.name) / "tempfile"
    try:
        yield tmp_fn
        tmp_fn.rename(output_path)
        tmp.cleanup()
    except:
        tmp.cleanup()
        raise


@contextlib.contextmanager
def intermediate_dir_with_cd(original_dir: Path, base_dir: Path = Path(".")):
    """
    Context manager to provide and change into a tempdir and ensure its removal and return to the
    original_dir upon exceptions.
    """

    if not original_dir.is_dir():
        raise DownloadError(f"Unable to setup temporary dir, original_dir does not exist: {original_dir}")

    tmp = tempfile.TemporaryDirectory(dir=base_dir)
    try:
        os.chdir(tmp.name)
        yield tmp
    finally:
        os.chdir(original_dir)
