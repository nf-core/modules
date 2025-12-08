import errno
import os

from refgenconf import MissingRecipeError
from ubiquerg import is_writable

from .asset_build_packages import asset_build_packages
from .exceptions import MissingFolderError


def _parse_user_build_input(input):
    """
    Parse user input specification. Used in build for specific parents and input parsing.

    :param Iterable[Iterable[str], ...] input: user command line input,
        formatted as follows: [[fasta=txt, test=txt], ...]
    :return dict: mapping of keys, which are input names and values
    """
    lst = []
    for i in input or []:
        lst.extend(i)
    return (
        {x.split("=")[0]: x.split("=")[1] for x in lst if "=" in x}
        if lst is not None
        else lst
    )


def _single_folder_writeable(d):
    return os.access(d, os.W_OK) and os.access(d, os.X_OK)


def _writeable(outdir, strict_exists=False):
    outdir = outdir or "."
    if os.path.exists(outdir):
        return _single_folder_writeable(outdir)
    elif strict_exists:
        raise MissingFolderError(outdir)
    return _writeable(os.path.dirname(outdir), strict_exists)


def _raise_missing_recipe_error(recipe):
    """
    Raise an error for a missing recipe, when one is requested

    :param str recipe: recipe name
    :raise MissingRecipeError: always
    """
    raise MissingRecipeError(
        f"Recipe '{recipe}' not found. Available recipes: "
        f"{', '.join(list(asset_build_packages.keys()))}"
    )


def _skip_lock(skip_arg, cfg):
    """
    If config read lock skip was not forced, check if dir is writable and set
    the default to the result

    :param bool skip_arg: argument selected on the CLI
    :param str cfg: path to the confjg
    :return bool: decision -- whether to skip the file lock for read
    """
    return is_writable(os.path.dirname(cfg)) if not skip_arg else True


def make_sure_path_exists(path):
    """
    Creates all directories in a path if it does not exist.

    :param str path: Path to create.
    :raises Exception: if the path creation attempt hits an error with
        a code indicating a cause other than pre-existence.
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
