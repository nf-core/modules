""" Helper functions """

import json
import logging
import os
import shutil
from copy import copy
from functools import partial
from re import sub
from typing import Iterable

from requests import ConnectionError, get
from ubiquerg import is_command_callable
from yacman import select_config

from .const import *
from .exceptions import DownloadJsonError, MissingAssetError
from .seqcol import SeqColClient

_LOGGER = logging.getLogger(__name__)

__all__ = ["select_genome_config", "get_dir_digest", "block_iter_repr"]


def select_genome_config(filename=None, conf_env_vars=CFG_ENV_VARS, **kwargs):
    """
    Get path to genome configuration file.

    :param str filename: name/path of genome configuration file
    :param Iterable[str] conf_env_vars: names of environment variables to
        consider; basically, a prioritized search list
    :return str: path to genome configuration file
    """
    return select_config(filename, conf_env_vars, **kwargs)


def unbound_env_vars(path):
    """
    Return collection of path parts that appear to be unbound env. vars.

    The given path is split on the active operating system's path delimiter;
    each resulting chunk is then deemed env.-var.-like if it begins with a
    dollar sign, and then os.getenv is queried to determine if it's bound.

    :param str path: Path to examine for unbound environment variables
    :return Iterable[str]: collection of path parts that appear to be unbound env. vars.
    """
    parts = path.split(os.path.sep)
    return [p for p in parts if p.startswith("$") and not os.getenv(p)]


def asciify_json_dict(json_dict):
    from ubiquerg.collection import asciify_dict

    return asciify_dict(json_dict)


def get_dir_digest(path, pm=None):
    """
    Generate a MD5 digest that reflects just the contents of the
    files in the selected directory.

    :param str path: path to the directory to digest
    :param pypiper.PipelineManager pm: a pipeline object, optional.
    The subprocess module will be used if not provided
    :return str: a digest, e.g. a3c46f201a3ce7831d85cf4a125aa334
    """
    if not is_command_callable("md5sum"):
        raise OSError(
            "md5sum command line tool is required for asset digest "
            "calculation. \n"
            "Install and try again, e.g on macOS: 'brew install "
            "md5sha1sum'"
        )
    cmd = (
        "cd {}; find . -type f -not -path './"
        + BUILD_STATS_DIR
        + "*' -exec md5sum {{}} \; | sort -k 2 | awk '{{print $1}}' | md5sum"
    )
    try:
        x = pm.checkprint(cmd.format(path))
    except AttributeError:
        try:
            from subprocess import check_output

            x = check_output(cmd.format(path), shell=True).decode("utf-8")
        except Exception as e:
            _LOGGER.warning(
                "{}: could not calculate digest for '{}'".format(
                    e.__class__.__name__, path
                )
            )
            return
    return str(sub(r"\W+", "", x))  # strips non-alphanumeric


def format_config_03_04(rgc, get_json_url):
    """
    upgrade the v0.3 config file format to v0.4 format:
    get genome digests from the server or local fasta assets,
    use the genome digests as primary key,
    add 'aliases' section to the config,
    remove 'genome_digests' section from the config
    replace all aliases in keys/asset names with genome digests

    :param obj rgc: RefGenConfV03 obj
    :param function(str, str) -> str get_json_url: how to build URL from
            genome server URL base, genome, and asset
    """

    _LOGGER.info("Upgrading v0.3 config file format to v0.4.")

    for genome, genome_v in rgc[CFG_GENOMES_KEY].items():
        digest = ""
        try:
            _LOGGER.info(
                f"Generating the digest from a local fasta file, "
                f"and creating the ASDs for {genome}."
            )
            tag = rgc.get_default_tag(genome, "fasta")
            asset_path = rgc.seek(genome, "fasta", tag, "fasta")
            ssc = SeqColClient({})
            digest, asdl = ssc.load_fasta(asset_path)
            _LOGGER.info(f"Generated {genome} digest from local fasta file: {digest}")
            # retrieve annotated sequence digests list to save in a JSON file
            pth = os.path.join(rgc[CFG_FOLDER_KEY], genome, genome + "__ASDs.json")
            os.makedirs(os.path.dirname(pth), exist_ok=True)
            with open(pth, "w") as jfp:
                json.dump(asdl, jfp)
            _LOGGER.info(f"Saved ASDs to JSON: {pth}")
        except (MissingAssetError, FileNotFoundError):
            _LOGGER.info(
                f"No local fasta asset found for {genome}. Retrieving digest from the server."
            )
            # get genome digest from the server
            cnt = 0
            servers = rgc[CFG_SERVERS_KEY]
            for server in servers:
                cnt += 1
                if not digest:
                    try:
                        url_alias = get_json_url(
                            s=server, i=API_VERSION + API_ID_ALIAS_DIGEST
                        ).format(alias=genome)
                        digest = send_data_request(url_alias)
                        _LOGGER.info(
                            f"Retrieved {genome} digest from the server: {digest}"
                        )
                    except (KeyError, ConnectionError, DownloadJsonError) as e:
                        if cnt == len(servers):
                            _LOGGER.info(
                                f"Failed to retrieve the digest for {genome}. "
                            )
                        continue
                continue

        if digest:
            # convert seek keys, children/parent asset keys from aliases to
            # genome digests
            rgc[CFG_GENOMES_KEY][genome] = replace_str_in_obj(genome_v, genome, digest)
            # use the genome digest as primary keys
            rgc[CFG_GENOMES_KEY][digest] = rgc[CFG_GENOMES_KEY].pop(genome)
            # create "aliases" section
            rgc[CFG_GENOMES_KEY][digest][CFG_ALIASES_KEY] = [genome]
            # remove old "genome_digest" section
            del rgc[CFG_GENOMES_KEY][digest][CFG_CHECKSUM_KEY]
        else:
            del rgc[CFG_GENOMES_KEY][genome]


def alter_file_tree_03_04(rgc, link_fun):
    """
    update file structure inside genome_folder:
    Drop genomes for which genome_digest is not available
    on any of the servers and do not have a fasta asset locally.
    contents inside genome_folder will be replaced by 'alias' and 'data' dir

    :param obj rgc: RefGenConfV03 obj
    :param callable link_fun: function to use to link files, e.g os.symlink
        or os.link
    """
    my_genome = {}
    for k, v in rgc[CFG_GENOMES_KEY].items():
        my_genome.update([(v[CFG_ALIASES_KEY][0], k)])

    _LOGGER.info(
        f"Creating '{DATA_DIR}' and '{ALIAS_DIR}' directories in "
        f"'{rgc[CFG_FOLDER_KEY]}'."
    )
    os.mkdir(os.path.abspath(os.path.join(rgc[CFG_FOLDER_KEY], DATA_DIR)))
    os.mkdir(os.path.abspath(os.path.join(rgc[CFG_FOLDER_KEY], ALIAS_DIR)))

    _LOGGER.info(
        f"Copying assets to '{DATA_DIR}' and creating alias symlinks in "
        f"'{ALIAS_DIR}'. Genomes that the digest could not be determined for "
        f"will be ignored."
    )
    for root, dirs, files in os.walk(rgc[CFG_FOLDER_KEY]):
        for dir in dirs:
            if dir in my_genome:
                shutil.copytree(
                    os.path.join(rgc[CFG_FOLDER_KEY], dir),
                    os.path.join(rgc[CFG_FOLDER_KEY], DATA_DIR, dir),
                    symlinks=True,
                )
        del dirs[:]

    for root, dirs, files in os.walk(os.path.join(rgc[CFG_FOLDER_KEY], DATA_DIR)):
        for dir in dirs:
            swap_names_in_tree(os.path.join(root, dir), my_genome[dir], dir)
            os.mkdir(os.path.join(rgc[CFG_FOLDER_KEY], ALIAS_DIR, dir))
            # create symlink for alias folder
            for genome, assets, files in os.walk(os.path.join(root, my_genome[dir])):
                for asset in assets:
                    old_path = os.path.join(genome, asset)
                    new_path = old_path.replace(my_genome[dir], dir).replace(
                        DATA_DIR, ALIAS_DIR
                    )
                    os.mkdir(new_path)

                for file in files:
                    des_path = os.path.join(genome, file)  # current file
                    src_path = (
                        os.path.realpath(des_path)
                        .replace(
                            os.path.realpath(rgc[CFG_FOLDER_KEY]),
                            os.path.join(rgc[CFG_FOLDER_KEY], DATA_DIR),
                        )
                        .replace(dir, my_genome[dir])
                    )  # replace /genome_folder with /genome_folder/data
                    # replace alias in the file name with genome digest

                    if os.path.islink(des_path):  # if the current file is a link
                        os.remove(
                            des_path
                        )  # remove the link that would not work after deleting old genome assest
                        link_fun(src_path, des_path)  # create the link with correct src

                    old_path = os.path.join(
                        genome, file
                    )  # path of the file in data dir
                    new_path = old_path.replace(
                        my_genome[dir], dir
                    ).replace(  # path of the file in alias
                        DATA_DIR, ALIAS_DIR
                    )

                    rel_old_path = os.path.join(
                        os.path.relpath(
                            os.path.dirname(old_path), os.path.dirname(new_path)
                        ),
                        os.path.basename(old_path),
                    )
                    link_fun(rel_old_path, new_path)
        del dirs[:]

    _LOGGER.info(
        f"Removing genome assets that have been copied " f"to '{DATA_DIR}' directory."
    )
    for genome, genome_v in rgc[CFG_GENOMES_KEY].items():
        d = os.path.join(rgc[CFG_FOLDER_KEY], genome_v[CFG_ALIASES_KEY][0])
        shutil.rmtree(d)


def swap_names_in_tree(top, new_name, old_name):
    """
    Rename all files and directories within a directory tree and the
    directory itself

    :param str top: path to the top of the tree to be renamed
    :param str new_name: new name
    :param str old_name: old name
    :return bool: whether the renaming has been carried out
    """

    def _rename(x, rt):
        os.rename(os.path.join(rt, x), os.path.join(rt, x.replace(old_name, new_name)))

    if not os.path.isdir(top):
        return False
    for root, dirs, files in os.walk(top):
        for directory in dirs:
            _rename(directory, root)
        for file in files:
            _rename(file, root)
    if os.path.split(top)[1] == old_name:
        # rename the top of the tree only if it is named as old_name
        os.rename(top, os.path.join(os.path.join(top, os.pardir), new_name))
    return True


def send_data_request(url, params=None):
    """
    Safely connect to the provided API endpoint and download the returned data.

    :param str url: server API endpoint
    :param dict params: query parameters
    :return dict: served data
    """
    _LOGGER.debug(f"Downloading JSON data; querying URL: {url}")
    resp = get(url, params=params)
    if resp.ok:
        try:
            return resp.json()
        except (json.JSONDecodeError, ValueError):
            _LOGGER.debug("The returned data is not a valid JSON")
            if resp.encoding == "utf-8" or resp.apparent_encoding == "ascii":
                _LOGGER.debug(f"Request returned pain text data: {resp.text}")
                return resp.text
    raise DownloadJsonError(resp)


def replace_str_in_obj(object, x, y):
    """
    Replace strings in an object

    :param any object: object to replace strings in
    :param str x: string to replace
    :param str y: replacement
    :return any: object with strings replaced
    """
    _replace = partial(replace_str_in_obj, x=x, y=y)
    obj = copy(object)
    if isinstance(obj, dict):
        for k, v in obj.items():
            obj[k] = _replace(v)
    if isinstance(obj, list):
        obj = [_replace(i) for i in obj]
    if isinstance(object, str):
        obj = object.replace(x, y)
    return obj


def block_iter_repr(input_obj, numbered=False):
    """
    Create a human readable string representation of an iterable. Either as a bulleted or numbered list.

    :param Iterable input_obj: object to create a representation for
    :param bool numbered: whether a numbered list should be created
    :param str: the representation
    """
    if isinstance(input_obj, str):
        input_obj = [input_obj]
    if not isinstance(input_obj, Iterable):
        raise TypeError("Input object has to be an Iterable")
    return (
        "\n"
        + "\n".join([" {}. {}".format(i + 1, val) for i, val in (enumerate(input_obj))])
        if numbered
        else "\n - {}".format("\n - ".join(input_obj))
    )
