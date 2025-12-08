#!/usr/bin/env python

import itertools
import json
import logging
import os
import re
import shutil
import signal
import sys
import warnings
from collections import OrderedDict
from collections.abc import Iterable, Mapping
from functools import partial
from inspect import getfullargspec as finspect
from urllib.error import ContentTooShortError, HTTPError
from urllib.parse import urlencode
from urllib.request import urlopen, urlretrieve

import yacman
from attmap import AttMap
from attmap import PathExAttMap as PXAM
from jsonschema.exceptions import ValidationError
from pkg_resources import iter_entry_points
from requests import ConnectionError
from requests.exceptions import MissingSchema
from rich.progress import BarColumn, Progress, TextColumn
from rich.table import Table
from ubiquerg import checksum, is_url, is_writable
from ubiquerg import parse_registry_path as prp
from ubiquerg import query_yes_no, untar

from .const import *
from .exceptions import *
from .helpers import (
    asciify_json_dict,
    block_iter_repr,
    get_dir_digest,
    select_genome_config,
    send_data_request,
)
from .progress_bar import _DownloadColumn, _TimeRemainingColumn, _TransferSpeedColumn
from .seqcol import SeqColClient

_LOGGER = logging.getLogger(__name__)

__all__ = ["RefGenConf", "upgrade_config"]


def _handle_sigint(filepath):
    def handle(sig, frame):
        _LOGGER.warning("\nThe download was interrupted: {}".format(filepath))
        try:
            os.remove(filepath)
        except OSError:
            _LOGGER.debug("'{}' not found, can't remove".format(filepath))
        else:
            _LOGGER.info("Incomplete file '{}' was removed".format(filepath))
        sys.exit(0)

    return handle


class RefGenConf(yacman.YacAttMap):
    """A sort of oracle of available reference genome assembly assets"""

    def __init__(
        self,
        filepath=None,
        entries=None,
        writable=False,
        wait_max=60,
        skip_read_lock=False,
        genome_exact=False,
        schema_source=None,
    ):
        """
        Create the config instance by with a filepath or key-value pairs.

        :param str filepath: a path to the YAML file to read
        :param Iterable[(str, object)] | Mapping[str, object] entries:
            config filepath or collection of key-value pairs
        :param bool writable: whether to create the object with write capabilities
        :param int wait_max: how long to wait for creating an object when the
            file that data will be read from is locked
        :param bool skip_read_lock: whether the file should not be locked for
            reading when object is created in read only mode
        :raise refgenconf.MissingConfigDataError: if a required configuration
            item is missing
        :raise ValueError: if entries is given as a string and is not a file
        """

        def _missing_key_msg(key, value):
            _LOGGER.debug("Config lacks '{}' key. Setting to: {}".format(key, value))

        super(RefGenConf, self).__init__(
            filepath=filepath,
            entries=entries,
            writable=writable,
            wait_max=wait_max,
            skip_read_lock=skip_read_lock,
            schema_source=schema_source or DEFAULT_CONFIG_SCHEMA,
            write_validate=True,
        )
        # assert correct config version
        try:
            version = self[CFG_VERSION_KEY]
        except KeyError:
            _missing_key_msg(CFG_VERSION_KEY, REQ_CFG_VERSION)
            self[CFG_VERSION_KEY] = REQ_CFG_VERSION
        else:
            try:
                version = float(version)
            except ValueError:
                _LOGGER.warning(
                    "Cannot parse config version as numeric: {}".format(version)
                )
            else:
                if version < REQ_CFG_VERSION:
                    msg = (
                        "This genome config (v{}) is not compliant with v{} standards. \n"
                        "To use current refgenconf, please use upgrade_config function to upgrade, or"
                        "downgrade refgenconf: 'pip install \"refgenconf>={},<{}\"'. \n"
                        "If refgenie is installed, you can use 'refgenie upgrade --target-version {}'".format(
                            self[CFG_VERSION_KEY],
                            str(REQ_CFG_VERSION),
                            REFGENIE_BY_CFG[str(version)],
                            REFGENIE_BY_CFG[str(REQ_CFG_VERSION)],
                            str(REQ_CFG_VERSION),
                        )
                    )
                    raise ConfigNotCompliantError(msg)

                else:
                    _LOGGER.debug("Config version is compliant: {}".format(version))

        # initialize "genomes_folder"
        if CFG_FOLDER_KEY not in self:
            self[CFG_FOLDER_KEY] = (
                os.path.dirname(filepath) if filepath else os.getcwd()
            )
            _missing_key_msg(CFG_FOLDER_KEY, self[CFG_FOLDER_KEY])
        # initialize "genome_servers"
        if CFG_SERVERS_KEY not in self and CFG_SERVER_KEY in self:
            # backwards compatibility after server config key change
            self[CFG_SERVERS_KEY] = self[CFG_SERVER_KEY]
            del self[CFG_SERVER_KEY]
            _LOGGER.debug(
                f"Moved servers list from '{CFG_SERVER_KEY}' to '{CFG_SERVERS_KEY}'"
            )
        try:
            if isinstance(self[CFG_SERVERS_KEY], list):
                tmp_list = [
                    server_url.rstrip("/") for server_url in self[CFG_SERVERS_KEY]
                ]
                self[CFG_SERVERS_KEY] = tmp_list
            else:  # Logic in pull_asset expects a list, even for a single server
                self[CFG_SERVERS_KEY] = self[CFG_SERVERS_KEY].rstrip("/")
                self[CFG_SERVERS_KEY] = [self[CFG_SERVERS_KEY]]
        except KeyError:
            _missing_key_msg(CFG_SERVERS_KEY, str([DEFAULT_SERVER]))
            self[CFG_SERVERS_KEY] = [DEFAULT_SERVER]

        # initialize "genomes" mapping
        if CFG_GENOMES_KEY in self:
            if not isinstance(self[CFG_GENOMES_KEY], PXAM):
                if self[CFG_GENOMES_KEY]:
                    _LOGGER.warning(
                        "'{k}' value is a {t_old}, not a {t_new}; setting to empty {t_new}".format(
                            k=CFG_GENOMES_KEY,
                            t_old=type(self[CFG_GENOMES_KEY]).__name__,
                            t_new=PXAM.__name__,
                        )
                    )
                self[CFG_GENOMES_KEY] = PXAM()
        else:
            self[CFG_GENOMES_KEY] = PXAM()

        self[CFG_GENOMES_KEY] = yacman.AliasedYacAttMap(
            entries=self[CFG_GENOMES_KEY],
            aliases=lambda x: {k: v.__getitem__(CFG_ALIASES_KEY) for k, v in x.items()},
            aliases_strict=True,
            exact=genome_exact,
        )

    def __bool__(self):
        minkeys = set(self.keys()) == set(RGC_REQ_KEYS)
        return not minkeys or bool(self[CFG_GENOMES_KEY])

    __nonzero__ = __bool__

    @property
    def plugins(self):
        """
        Plugins registered by entry points in the current Python env

        :return dict[dict[function(refgenconf.RefGenConf)]]: dict which keys
            are names of all possible hooks and values are dicts mapping
            registered functions names to their values
        """
        return {
            h: {ep.name: ep.load() for ep in iter_entry_points("refgenie.hooks." + h)}
            for h in HOOKS
        }

    @property
    def genome_aliases(self):
        """
        Mapping of human-readable genome identifiers to genome identifiers

        :return dict: mapping of human-readable genome identifiers to genome
            identifiers
        """
        return self.genomes[yacman.IK][yacman.ALIASES_KEY]

    @property
    def genome_aliases_table(self):
        """
        Mapping of human-readable genome identifiers to genome identifiers

        :return dict: mapping of human-readable genome identifiers to genome
            identifiers
        """
        table = Table(title="Genome aliases")
        table.add_column("genome")
        table.add_column("alias")
        if CFG_GENOMES_KEY not in self or not self[CFG_GENOMES_KEY]:
            return table
        for genome, genome_dict in self[CFG_GENOMES_KEY].items():
            if (
                CFG_ALIASES_KEY not in self[CFG_GENOMES_KEY][genome]
                or not self[CFG_GENOMES_KEY][genome][CFG_ALIASES_KEY]
            ):
                aliases = ""
            else:
                aliases = ", ".join(self[CFG_GENOMES_KEY][genome][CFG_ALIASES_KEY])
            table.add_row(genome, aliases)
        return table

    @property
    def data_dir(self):
        """
        Path to the genome data directory

        :return str: path to the directory where the assets are stored
        """
        return os.path.abspath(os.path.join(self[CFG_FOLDER_KEY], DATA_DIR))

    @property
    def alias_dir(self):
        """
        Path to the genome alias directory

        :return str: path to the directory where the assets are stored
        """
        return os.path.abspath(os.path.join(self[CFG_FOLDER_KEY], ALIAS_DIR))

    @property
    def file_path(self):
        """
        Path to the genome configuration file

        :return str: path to the genome configuration file
        """
        try:
            return self[yacman.IK][yacman.FILEPATH_KEY]
        except (AttributeError, KeyError):
            return None

    def initialize_config_file(self, filepath=None):
        """
        Initialize genome configuration file on disk

        :param str filepath: a valid path where the configuration file should be initialized
        :return str: the filepath the file was initialized at
        :raise OSError: in case the file could not be initialized due to insufficient permissions or pre-existence
        :raise TypeError: if no valid filepath cat be determined
        """

        def _write_fail_err(reason):
            raise OSError("Can't initialize, {}: {} ".format(reason, filepath))

        filepath = select_genome_config(filepath, check_exist=False)
        if not isinstance(filepath, str):
            raise TypeError(
                f"Could not determine a valid path to initialize a "
                f"configuration file: {filepath}"
            )
        if os.path.exists(filepath):
            _write_fail_err("file exists")
        if not is_writable(filepath, check_exist=False):
            _write_fail_err("insufficient permissions")
        self.make_writable(filepath)
        self.write()
        self.make_readonly()
        _LOGGER.info(f"Initialized genome configuration file: {filepath}")
        os.makedirs(self.data_dir, exist_ok=True)
        os.makedirs(self.alias_dir, exist_ok=True)
        _LOGGER.info(
            f"Created directories:{block_iter_repr([self.data_dir, self.alias_dir])}"
        )

        return filepath

    def list(self, genome=None, order=None, include_tags=False):
        """
        List local assets; map each namespace to a list of available asset names

        :param callable(str) -> object order: how to key genome IDs for sort
        :param list[str] | str genome: genomes that the assets should be found for
        :param bool include_tags: whether asset tags should be included in the returned dict
        :return Mapping[str, Iterable[str]]: mapping from assembly name to
            collection of available asset names.
        """
        self.run_plugins(PRE_LIST_HOOK)
        refgens = self._select_genomes(genome=genome, order=order)
        if include_tags:
            self.run_plugins(POST_LIST_HOOK)
            return OrderedDict(
                [
                    (
                        g,
                        sorted(
                            _make_asset_tags_product(
                                self[CFG_GENOMES_KEY][g][CFG_ASSETS_KEY], ":"
                            ),
                            key=order,
                        ),
                    )
                    for g in refgens
                    if CFG_ASSETS_KEY in self[CFG_GENOMES_KEY][g]
                ]
            )
        self.run_plugins(POST_LIST_HOOK)
        return OrderedDict(
            [
                (
                    g,
                    sorted(
                        list(self[CFG_GENOMES_KEY][g][CFG_ASSETS_KEY].keys()), key=order
                    ),
                )
                for g in refgens
                if CFG_ASSETS_KEY in self[CFG_GENOMES_KEY][g]
            ]
        )

    def get_asset_table(
        self,
        genomes=None,
        server_url=None,
        get_json_url=lambda s, i: construct_request_url(s, i, PRIVATE_API),
    ):
        """
        Get a rich.Table object representing assets available locally

        :param list[str] genomes: genomes to restrict the results with
        :param str server_url: server URL to query for the remote genome data
        :param function(str, str) -> str get_json_url: how to build URL from
            genome server URL base, genome, and asset
        :return rich.table.Table: table of assets available locally
        """

        def _fill_table_with_genomes_data(rgc, genomes_data, table, genomes=None):
            it = "([italic]{}[/italic])"
            table.add_column("genome")
            if genomes:
                table.add_column("asset " + it.format("seek_keys"))
                table.add_column("tags")
                for g in genomes:
                    try:
                        genome = rgc.get_genome_alias_digest(alias=g, fallback=True)
                    except yacman.UndefinedAliasError:
                        rgc.set_genome_alias(
                            genome=g, create_genome=True, no_write=True
                        )
                        genome = rgc.get_genome_alias_digest(alias=g, fallback=True)
                    if genome not in genomes_data:
                        _LOGGER.error(f"Genome {g} ({genome}) not found")
                        continue
                    genome_dict = genomes_data[genome]
                    if CFG_ASSETS_KEY not in genome_dict:
                        continue
                    for asset, asset_dict in genome_dict[CFG_ASSETS_KEY].items():
                        tags = list(asset_dict[CFG_ASSET_TAGS_KEY].keys())
                        if (
                            CFG_SEEK_KEYS_KEY
                            not in asset_dict[CFG_ASSET_TAGS_KEY][tags[0]]
                        ):
                            continue
                        seek_keys = list(
                            asset_dict[CFG_ASSET_TAGS_KEY][tags[0]][
                                CFG_SEEK_KEYS_KEY
                            ].keys()
                        )
                        table.add_row(
                            ", ".join(genome_dict[CFG_ALIASES_KEY]),
                            "{} ".format(asset) + it.format(", ".join(seek_keys)),
                            ", ".join(tags),
                        )
            else:
                table.add_column("assets")
                for genome in list(genomes_data.keys()):
                    genome_dict = genomes_data[genome]
                    if CFG_ASSETS_KEY not in genome_dict:
                        continue
                    table.add_row(
                        ", ".join(genome_dict[CFG_ALIASES_KEY]),
                        ", ".join(list(genome_dict[CFG_ASSETS_KEY].keys())),
                    )
            return table

        if server_url is None:
            genomes_data = self[CFG_GENOMES_KEY]
            title = (
                f"Local refgenie assets\nServer subscriptions: "
                f"{', '.join(self[CFG_SERVERS_KEY])}"
            )
        else:
            genomes_data = send_data_request(
                get_json_url(server_url, API_ID_GENOMES_DICT)
            )
            title = f"Remote refgenie assets\nServer URL: {server_url}"
        c = (
            f"use refgenie list{'r' if server_url is not None else ''} "
            f"-g <genome> for more detailed view"
            if genomes is None
            else ""
        )
        return _fill_table_with_genomes_data(
            self, genomes_data, Table(title=title, min_width=70, caption=c), genomes
        )

    def assets_str(
        self,
        offset_text="  ",
        asset_sep=", ",
        genome_assets_delim="/ ",
        genome=None,
        order=None,
    ):
        """
        Create a block of text representing genome-to-asset mapping.

        :param str offset_text: text that begins each line of the text
            representation that's produced
        :param str asset_sep: the delimiter between names of types of assets,
            within each genome line
        :param str genome_assets_delim: the delimiter to place between
            reference genome assembly name and its list of asset names
        :param list[str] | str genome: genomes that the assets should be found for
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return str: text representing genome-to-asset mapping
        """
        refgens = self._select_genomes(genome=genome, order=order)
        make_line = partial(
            _make_genome_assets_line,
            offset_text=offset_text,
            genome_assets_delim=genome_assets_delim,
            asset_sep=asset_sep,
            order=order,
            rjust=max(map(len, refgens) or [0]) + 2,
        )
        return "\n".join(
            [make_line(g, self[CFG_GENOMES_KEY][g][CFG_ASSETS_KEY]) for g in refgens]
        )

    def add(self, path, genome, asset, tag=None, seek_keys=None, force=False):
        """
        Add an external asset to the config

        :param str path: a path to the asset to add; must exist and be relative
            to the genome_folder
        :param str genome: genome name
        :param str asset: asset name
        :param str tag: tag name
        :param dict seek_keys: seek keys to add
        :param bool force: whether to force existing asset overwrite
        """
        try:
            genome = self.get_genome_alias_digest(alias=genome, fallback=True)
        except yacman.UndefinedAliasError:
            _LOGGER.error(
                "No digest defined for '{}'. Set an alias or pull an"
                " asset to initialize.".format(genome)
            )
            return False
        tag = tag or self.get_default_tag(genome, asset)
        abspath = os.path.join(self[CFG_FOLDER_KEY], path)
        remove = False
        if not os.path.exists(abspath) or not os.path.isabs(abspath):
            raise OSError(
                "Provided path must exist and be relative to the"
                " genome_folder: {}".format(self[CFG_FOLDER_KEY])
            )
        try:
            _assert_gat_exists(self[CFG_GENOMES_KEY], genome, asset, tag)
        except Exception:
            pass
        else:
            if not force and not query_yes_no(
                "'{}/{}:{}' exists. Do you want to overwrite?".format(
                    genome, asset, tag
                )
            ):
                _LOGGER.info("Aborted by a user, asset no added")
                return False
            remove = True
            _LOGGER.info("Will remove existing to overwrite")
        tag_data = {
            CFG_ASSET_PATH_KEY: path,
            CFG_ASSET_CHECKSUM_KEY: get_dir_digest(abspath) or "",
        }
        msg = "Added asset: {}/{}:{} {}".format(
            genome,
            asset,
            tag,
            "" if not seek_keys else "with seek keys: {}".format(seek_keys),
        )
        if not self.file_path:
            if remove:
                self.cfg_remove_assets(genome, asset, tag)
            self.update_tags(genome, asset, tag, tag_data)
            self.update_seek_keys(genome, asset, tag, seek_keys or {asset: "."})
            self.set_default_pointer(genome, asset, tag)
            _LOGGER.info(msg)
        else:
            with self as rgc:
                if remove:
                    rgc.cfg_remove_assets(genome, asset, tag)
                rgc.update_tags(genome, asset, tag, tag_data)
                rgc.update_seek_keys(genome, asset, tag, seek_keys or {asset: "."})
                rgc.set_default_pointer(genome, asset, tag)
                _LOGGER.info(msg)
        self._symlink_alias(genome, asset, tag)
        return True

    def get_symlink_paths(self, genome, asset=None, tag=None, all_aliases=False):
        """
        Get path to the alias directory for the selected genome-asset-tag

        :param str genome: reference genome ID
        :param str asset: asset name
        :param str tag: tag name
        :param bool all_aliases: whether to return a collection of symbolic
            links or just the first one from the alias list
        :return dict:
        """
        try:
            defined_aliases = self.get_genome_alias(
                genome, fallback=True, all_aliases=all_aliases
            )
        except yacman.UndefinedAliasError:
            return {}
        alias = _make_list_of_str(defined_aliases)
        if asset:
            tag = tag or self.get_default_tag(genome, asset)
        return {
            a: os.path.join(self.alias_dir, a, asset, tag)
            if asset
            else os.path.join(self.alias_dir, a)
            for a in alias
        }

    def _symlink_alias(
        self, genome, asset=None, tag=None, link_fun=lambda t, s: os.symlink(t, s)
    ):
        """
        Go through the files in the asset directory and recreate the asset
        directory tree, but instead of copying files, create symbolic links

        :param str genome: reference genome ID
        :param str asset: asset name
        :param str tag: tag name
        :param callable link_fun: function to use to link files, e.g os.symlink
            or os.link
        """

        def _rpl(x):
            """
            RePLace genome digest with human-readable genome ID, if exists

            :param str x: string to replace digest with alias in
            :return str: processed string
            """
            return x.replace(genome_digest, alias)

        if not callable(link_fun) or len(finspect(link_fun).args) != 2:
            raise TypeError(
                "Linking function must be a two-arg function (target, destination)"
            )
        created = []
        genome_digest = self.get_genome_alias_digest(genome, fallback=True)
        if asset:
            tag = tag or self.get_default_tag(genome, asset)
            src_path = self.seek_src(genome, asset, tag, enclosing_dir=True)
        else:
            src_path = os.path.join(self.data_dir, genome_digest)
        target_paths_mapping = self.get_symlink_paths(
            genome_digest, asset, tag, all_aliases=True
        )
        for alias, path in target_paths_mapping.items():
            if not os.path.exists(path):
                os.makedirs(path, exist_ok=True)
                for root, dirs, files in os.walk(src_path):
                    appendix = os.path.relpath(root, src_path)
                    for directory in dirs:
                        try:
                            os.makedirs(os.path.join(path, appendix, _rpl(directory)))
                        except FileExistsError:
                            continue
                    for file in files:
                        try:
                            rel = os.path.relpath(
                                os.path.join(root, file), os.path.join(path, appendix)
                            )
                            new_path = os.path.join(path, appendix, _rpl(file))
                            link_fun(rel, new_path)
                        except FileExistsError:
                            _LOGGER.warning(
                                f"Could not create link, file exists: {new_path}"
                            )
                            continue
                created.append(path)
        if created:
            _LOGGER.info(f"Created alias directories:{block_iter_repr(created)}")

    @staticmethod
    def _remove_symlink_alias(symlink_dict, aliases_to_remove):
        """
        Remove the symlink directories

        :param list[str] | str aliases_to_remove: collection of aliases to
            remove the symlink directories for
        :param dict symlink_dict: a dictionary mapping alias names to the
            respective symlink directories
        """
        dirs_to_remove = [symlink_dict[k] for k in _make_list_of_str(aliases_to_remove)]
        for d in dirs_to_remove:
            shutil.rmtree(d)
        if dirs_to_remove:
            _LOGGER.info(
                "Removed alias directories: \n - {}".format(
                    "\n - ".join(dirs_to_remove)
                )
            )

    def filepath(self, genome, asset, tag, ext=".tgz", dir=False):
        """
        Determine path to a particular asset for a particular genome.

        :param str genome: reference genome ID
        :param str asset: asset name
        :param str tag: tag name
        :param str ext: file extension
        :param bool dir: whether to return the enclosing directory instead of the file
        :return str: path to asset for given genome and asset kind/name
        """
        tag_dir = os.path.join(self.data_dir, genome, asset, tag)
        return os.path.join(tag_dir, asset + "__" + tag + ext) if not dir else tag_dir

    def genomes_list(self, order=None):
        """
        Get a list of this configuration's reference genome assembly IDs.

        :return Iterable[str]: list of this configuration's reference genome
            assembly IDs
        """
        return sorted(
            [
                self.get_genome_alias(x, fallback=True)
                for x in self[CFG_GENOMES_KEY].keys()
            ],
            key=order,
        )

    def genomes_str(self, order=None):
        """
        Get as single string this configuration's reference genome assembly IDs.

        :param function(str) -> object order: how to key genome IDs for sort
        :return str: single string that lists this configuration's known
            reference genome assembly IDs
        """
        return ", ".join(self.genomes_list(order))

    def seek(
        self,
        genome_name,
        asset_name,
        tag_name=None,
        seek_key=None,
        strict_exists=None,
        enclosing_dir=False,
        all_aliases=False,
        check_exist=lambda p: os.path.exists(p) or is_url(p),
    ):
        """
        Seek path to a specified genome-asset-tag alias

        :param str genome_name: name of a reference genome assembly of interest
        :param str asset_name: name of the particular asset to fetch
        :param str tag_name: name of the particular asset tag to fetch
        :param str seek_key: name of the particular subasset to fetch
        :param bool | NoneType strict_exists: how to handle case in which
            path doesn't exist; True to raise IOError, False to raise
            RuntimeWarning, and None to do nothing at all.
            Default: None (do not check).
        :param function(callable) -> bool check_exist: how to check for
            asset/path existence
        :param bool enclosing_dir: whether a path to the entire enclosing
            directory should be returned, e.g. for a fasta asset that has 3
            seek_keys pointing to 3 files in an asset dir, that asset dir
            is returned
        :param bool all_aliases: whether to return paths to all asset aliases or
            just the one for the specified 'genome_name` argument
        :return str: path to the asset
        :raise TypeError: if the existence check is not a one-arg function
        :raise refgenconf.MissingGenomeError: if the named assembly isn't known
            to this configuration instance
        :raise refgenconf.MissingAssetError: if the names assembly is known to
            this configuration instance, but the requested asset is unknown
        """
        tag_name = tag_name or self.get_default_tag(genome_name, asset_name)
        try:
            genome_digest = self.get_genome_alias_digest(genome_name, fallback=True)
        except yacman.UndefinedAliasError:
            raise MissingGenomeError(f"Your genomes do not include '{genome_name}'")
        genome_ids = _make_list_of_str(
            self.get_genome_alias(genome_digest, fallback=True, all_aliases=True)
        )
        idx = 0
        if genome_name in genome_ids:
            idx = genome_ids.index(genome_name)
        self._assert_gat_exists(genome_name, asset_name, tag_name)
        asset_tag_data = self[CFG_GENOMES_KEY][genome_name][CFG_ASSETS_KEY][asset_name][
            CFG_ASSET_TAGS_KEY
        ][tag_name]
        if not seek_key:
            if asset_name in asset_tag_data[CFG_SEEK_KEYS_KEY]:
                seek_val = asset_tag_data[CFG_SEEK_KEYS_KEY][asset_name]
            else:
                seek_val = ""
        else:
            try:
                seek_val = asset_tag_data[CFG_SEEK_KEYS_KEY][seek_key]
            except KeyError:
                if seek_key == "dir":
                    seek_val = "."
                else:
                    raise MissingSeekKeyError(
                        f"Seek key '{seek_key}' not defined for: "
                        f"'{genome_name}.{asset_name}:{tag_name}'"
                    )
        if enclosing_dir:
            seek_val = ""
        fullpath = os.path.join(
            self.alias_dir, genome_digest, asset_name, tag_name, seek_val
        )
        fullpaths = [fullpath.replace(genome_digest, gid) for gid in genome_ids]
        paths_existence = [check_exist(fp) for fp in fullpaths]
        if all(paths_existence):
            return fullpaths if all_aliases else fullpaths[idx]
        nonexistent_pths = [
            fullpaths[p] for p in [i for i, x in enumerate(paths_existence) if not x]
        ]
        msg = "For genome '{}' path to the asset '{}/{}:{}' doesn't exist: {}".format(
            genome_name,
            genome_name,
            asset_name,
            seek_key,
            tag_name,
            ", ".join(nonexistent_pths),
        )
        if strict_exists is None:
            _LOGGER.debug(msg)
        elif strict_exists is True:
            raise OSError(msg)
        else:
            warnings.warn(msg, RuntimeWarning)
        return fullpaths if all_aliases else fullpaths[idx]

    def seekr(
        self,
        genome_name,
        asset_name,
        tag_name=None,
        seek_key=None,
        remote_class="http",
        get_url=lambda server, id: construct_request_url(server, id),
    ):
        """
        Seek a remote path to a specified genome/asset.seek_key:tag

        :param str genome_name: name of a reference genome assembly of interest
        :param str asset_name: name of the particular asset to fetch
        :param str tag_name: name of the particular asset tag to fetch
        :param str seek_key: name of the particular subasset to fetch
        :param str remote_class: remote data provider class, e.g. 'http' or 's3'
        :param function(serverUrl, operationId) -> str get_url: how to determine
            URL request, given server URL and endpoint operationID
        :return str: path to the asset
        """
        good_servers = [
            s for s in self[CFG_SERVERS_KEY] if get_url(s, API_ID_ASSET_PATH)
        ]
        _LOGGER.debug(f"Compatible refgenieserver instances: {good_servers}")
        for url in good_servers:
            try:
                genome_digest = self.get_genome_alias_digest(alias=genome_name)
            except yacman.UndefinedAliasError:
                _LOGGER.info(f"No local digest for genome alias: {genome_name}")
                if not self.set_genome_alias(
                    genome=genome_name, servers=[url], create_genome=True
                ):
                    continue
                genome_digest = self.get_genome_alias_digest(alias=genome_name)

            asset_seek_key_url = get_url(url, API_ID_ASSET_PATH).format(
                genome=genome_digest, asset=asset_name, seek_key=seek_key or asset_name
            )
            if asset_seek_key_url is None:
                continue
            asset_seek_key_target = send_data_request(
                asset_seek_key_url,
                params={"tag": tag_name, "remoteClass": remote_class},
            )
            return asset_seek_key_target

    def seek_src(
        self,
        genome_name,
        asset_name,
        tag_name=None,
        seek_key=None,
        strict_exists=None,
        enclosing_dir=False,
        check_exist=lambda p: os.path.exists(p) or is_url(p),
    ):
        """
        Seek path to a specified genome-asset-tag

        :param str genome_name: name of a reference genome assembly of interest
        :param str asset_name: name of the particular asset to fetch
        :param str tag_name: name of the particular asset tag to fetch
        :param str seek_key: name of the particular subasset to fetch
        :param bool | NoneType strict_exists: how to handle case in which
            path doesn't exist; True to raise IOError, False to raise
            RuntimeWarning, and None to do nothing at all.
            Default: None (do not check).
        :param function(callable) -> bool check_exist: how to check for
            asset/path existence
        :param bool enclosing_dir: whether a path to the entire enclosing
            directory should be returned, e.g. for a fasta asset that has 3
            seek_keys pointing to 3 files in an asset dir, that asset dir
            is returned
        :return str: path to the asset
        :raise TypeError: if the existence check is not a one-arg function
        :raise refgenconf.MissingGenomeError: if the named assembly isn't known
            to this configuration instance
        :raise refgenconf.MissingAssetError: if the names assembly is known to
            this configuration instance, but the requested asset is unknown
        """
        tag_name = tag_name or self.get_default_tag(genome_name, asset_name)
        _LOGGER.debug(
            "getting asset: '{}/{}.{}:{}'".format(
                genome_name, asset_name, seek_key, tag_name
            )
        )
        if not callable(check_exist) or len(finspect(check_exist).args) != 1:
            raise TypeError("Asset existence check must be a one-arg function.")
        # 3 'path' key options supported
        # option1: absolute path
        # get just the saute path value from the config
        path_val = _genome_asset_path(
            self[CFG_GENOMES_KEY],
            genome_name,
            asset_name,
            tag_name,
            enclosing_dir=True,
            no_tag=True,
            seek_key=None,
        )
        _LOGGER.debug("Trying absolute path: {}".format(path_val))
        if seek_key:
            path = os.path.join(path_val, seek_key)
        else:
            path = path_val
        if os.path.isabs(path) and check_exist(path):
            return path
        genome_name = self.get_genome_alias_digest(genome_name, fallback=True)
        # option2: relative to genome_folder/{genome} (default, canonical)
        path = _genome_asset_path(
            self[CFG_GENOMES_KEY],
            genome_name,
            asset_name,
            tag_name,
            seek_key,
            enclosing_dir,
        )
        fullpath = os.path.join(self.data_dir, genome_name, path)
        _LOGGER.debug(
            "Trying relative to genome_folder/genome/_data ({}/{}/{}): {}".format(
                self[CFG_FOLDER_KEY], genome_name, DATA_DIR, fullpath
            )
        )
        if check_exist(fullpath):
            return fullpath
        # option3: relative to the genome_folder (if option2 does not exist)
        gf_relpath = os.path.join(
            self[CFG_FOLDER_KEY],
            _genome_asset_path(
                self[CFG_GENOMES_KEY],
                genome_name,
                asset_name,
                tag_name,
                seek_key,
                enclosing_dir,
                no_tag=True,
            ),
        )
        _LOGGER.debug(
            "Trying path relative to genome_folder ({}): {}".format(
                self[CFG_FOLDER_KEY], gf_relpath
            )
        )
        if check_exist(gf_relpath):
            return gf_relpath

        msg = "For genome '{}' the asset '{}.{}:{}' doesn't exist; tried: {}".format(
            genome_name,
            asset_name,
            seek_key,
            tag_name,
            ", ".join([path, gf_relpath, fullpath]),
        )
        # return option2 if existence not enforced
        if strict_exists is None:
            _LOGGER.debug(msg)
        elif strict_exists is True:
            raise OSError(msg)
        else:
            warnings.warn(msg, RuntimeWarning)
        return fullpath

    def get_default_tag(self, genome, asset, use_existing=True):
        """
        Determine the asset tag to use as default. The one indicated by
        the 'default_tag' key in the asset section is returned.
        If no 'default_tag' key is found, by default the first listed tag is returned
        with a RuntimeWarning. This behavior can be turned off with use_existing=False

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of the particular asset of interest
        :param bool use_existing: whether the first tag in the config should be
            returned in case there is no default tag defined for an asset
        :return str: name of the tag to use as the default one
        """
        try:
            self._assert_gat_exists(genome, asset)
        except RefgenconfError:
            _LOGGER.info(
                "Using '{}' as the default tag for '{}/{}'".format(
                    DEFAULT_TAG, genome, asset
                )
            )
            return DEFAULT_TAG
        try:
            return self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                CFG_ASSET_DEFAULT_TAG_KEY
            ]
        except KeyError:
            alt = (
                self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                    CFG_ASSET_TAGS_KEY
                ].keys()[0]
                if use_existing
                else DEFAULT_TAG
            )
            if isinstance(alt, str):
                if alt != DEFAULT_TAG:
                    warnings.warn(
                        "Could not find the '{}' key for asset '{}/{}'. "
                        "Used the first one in the config instead: '{}'. "
                        "Make sure it does not corrupt your workflow.".format(
                            CFG_ASSET_DEFAULT_TAG_KEY, genome, asset, alt
                        ),
                        RuntimeWarning,
                    )
                else:
                    warnings.warn(
                        "Could not find the '{}' key for asset '{}/{}'. Returning '{}' "
                        "instead. Make sure it does not corrupt your workflow.".format(
                            CFG_ASSET_DEFAULT_TAG_KEY, genome, asset, alt
                        ),
                        RuntimeWarning,
                    )
                return alt
        except TypeError:
            _raise_not_mapping(
                self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset], "Asset section "
            )

    def set_default_pointer(
        self,
        genome,
        asset,
        tag,
        force_exists=False,
        force_digest=None,
        force_fasta=False,
    ):
        """
        Point to the selected tag by default

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of the particular asset of interest
        :param str tag: name of the particular asset tag to point to by default
        :param str force_digest: digest to force update of. The alias will
            not be converted to the digest, even if provided.
        :param bool force_fasta: whether setting a default tag for a fasta asset
            should be forced. Beware: This could lead to genome identity issues
        :param bool force_exists: whether the default tag change should be
            forced (even if it exists)
        """
        self._assert_gat_exists(genome, asset, tag)
        asset_dict = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset]
        if (
            CFG_ASSET_DEFAULT_TAG_KEY in asset_dict
            and len(asset_dict[CFG_ASSET_DEFAULT_TAG_KEY]) > 0
        ):
            if not force_exists:
                return
            if asset == "fasta" and not force_fasta:
                raise NotImplementedError(
                    "Can't change the default tag for fasta assets, "
                    "this would lead to genome identity issues"
                )
        self.update_assets(
            genome, asset, {CFG_ASSET_DEFAULT_TAG_KEY: tag}, force_digest=force_digest
        )
        _LOGGER.info(f"Default tag for '{genome}/{asset}' set to: {tag}")

    def list_assets_by_genome(self, genome=None, order=None, include_tags=False):
        """
        List types/names of assets that are available for one--or all--genomes.

        :param str | NoneType genome: reference genome assembly ID, optional;
            if omitted, the full mapping from genome to asset names
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :param bool include_tags: whether asset tags should be included in the
            returned dict
        :return Iterable[str] | Mapping[str, Iterable[str]]: collection of
            asset type names available for particular reference assembly if
            one is provided, else the full mapping between assembly ID and
            collection available asset type names
        """
        if genome:
            genome = self.get_genome_alias(digest=genome, fallback=True)
        return (
            self.list(genome, order, include_tags=include_tags)[genome]
            if genome is not None
            else self.list(order, include_tags=include_tags)
        )

    def list_genomes_by_asset(self, asset=None, order=None):
        """
        List assemblies for which a particular asset is available.

        :param str | NoneType asset: name of type of asset of interest, optional
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return Iterable[str] | Mapping[str, Iterable[str]]: collection of
            assemblies for which the given asset is available; if asset
            argument is omitted, the full mapping from name of asset type to
            collection of assembly names for which the asset key is available
            will be returned.
        """
        return (
            self._invert_genomes(order)
            if not asset
            else sorted(
                [
                    self.get_genome_alias(g, fallback=True)
                    for g, data in self[CFG_GENOMES_KEY].items()
                    if asset in data.get(CFG_ASSETS_KEY)
                ],
                key=order,
            )
        )

    def list_seek_keys_values(self, genomes=None, assets=None):
        """
        List values for all seek keys for the specified genome and asset.
        Leave the arguments out to get all seek keys values managed by refgenie.

        :param str | List[str] genome_names: optional list of genomes to include
        :param str | List[str] asset_names: optional list of assets to include
        :return dict: a nested dictionary with the seek key values
        """
        ret = {}

        if genomes is None:
            genome_names = self.genomes_list()
        else:
            genome_names = _make_list_of_str(genomes)

        for genome_name in genome_names:
            self._assert_gat_exists(genome_name)
            ret[genome_name] = {}
            if assets is None:
                asset_names = self.list_assets_by_genome(genome_name)
            else:
                asset_names = _make_list_of_str(assets)
            for asset_name in asset_names:
                try:
                    self._assert_gat_exists(genome_name, asset_name)
                except MissingAssetError as e:
                    _LOGGER.warning(f"Skipping {asset_name} asset: {str(e)}")
                    continue
                asset_mapping = self[CFG_GENOMES_KEY][genome_name][CFG_ASSETS_KEY][
                    asset_name
                ]
                ret[genome_name][asset_name] = {}
                for tag_name in get_asset_tags(asset_mapping):
                    tag_mapping = asset_mapping[CFG_ASSET_TAGS_KEY][tag_name]
                    ret[genome_name][asset_name][tag_name] = {}
                    for seek_key_name in get_tag_seek_keys(tag_mapping):
                        ret[genome_name][asset_name][tag_name][
                            seek_key_name
                        ] = self.seek(genome_name, asset_name, tag_name, seek_key_name)
        return ret

    def get_local_data_str(self, genome=None, order=None):
        """
        List locally available reference genome IDs and assets by ID.

        :param list[str] | str genome: genomes that the assets should be found for
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return str, str: text reps of locally available genomes and assets
        """
        exceptions = []
        if genome is not None:
            genome = _make_list_of_str(genome)
            for g in genome:
                try:
                    self._assert_gat_exists(gname=g)
                except MissingGenomeError as e:
                    exceptions.append(e)
            if exceptions:
                raise MissingGenomeError(", ".join(map(str, exceptions)))
        return (
            ", ".join(self._select_genomes(genome=genome, order=order)),
            self.assets_str(genome=genome, order=order),
        )

    def get_remote_data_str(
        self,
        genome=None,
        order=None,
        get_url=lambda server, id: construct_request_url(server, id),
    ):
        """
        List genomes and assets available remotely.

        :param function(serverUrl, operationId) -> str get_url: how to determine
            URL request, given server URL and endpoint operationID
        :param list[str] | str genome: genomes that the assets should be found for
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return str, str: text reps of remotely available genomes and assets
        """
        warnings.warn(
            "Please use listr method instead; get_remote_data_str will be "
            "removed in the next release.",
            category=DeprecationWarning,
        )
        return self.listr(genome, order, get_url)

    def listr(
        self,
        genome=None,
        get_url=lambda server, id: construct_request_url(server, id),
        as_digests=False,
    ):
        """
        List genomes and assets available remotely on all servers the object
        subscribes to

        :param function(serverUrl, operationId) -> str get_url: how to determine
            URL request, given server URL and endpoint operationID
        :param list[str] | str genome: genomes that the assets should be found for
        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return dict[OrderedDict[list]]: remotely available genomes and assets
            keyed by genome keyed by source server endpoint
        """
        data_by_server = {}

        for url in self[CFG_SERVERS_KEY]:
            aliases_url = get_url(url, API_ID_ALIASES_DICT)
            assets_url = get_url(url, API_ID_ASSETS)
            if assets_url is None or aliases_url is None:
                continue

            aliases_by_digest = send_data_request(aliases_url)
            # convert the original, condensed mapping to a data structure with optimal time complexity
            digests_by_alias = {}
            for k, v in aliases_by_digest.items():
                for alias in v:
                    digests_by_alias[alias] = k

            genome_digests = None
            genomes = genome if isinstance(genome, list) else [genome]
            if genome is not None:
                genome_digests = [
                    g
                    if g in aliases_by_digest.keys()
                    else digests_by_alias.get(g, None)
                    for g in genomes
                ]
                if genome_digests is None:
                    _LOGGER.info(f"{genome} not found on server: {url}")
                    continue

            server_data = self._list_remote(
                url=assets_url,
                genome=genome_digests,
            )
            data_by_server[assets_url] = (
                server_data
                if as_digests
                else {aliases_by_digest[k][0]: v for k, v in server_data.items()}
            )

        return data_by_server

    def tag(self, genome, asset, tag, new_tag, files=True, force=False):
        """
        Retags the asset selected by the tag with the new_tag.
        Prompts if default already exists and overrides upon confirmation.

        This method does not override the original asset entry in the RefGenConf
        object. It creates its copy and tags it with the new_tag.
        Additionally, if the retagged asset has any children their parent will
        be retagged as new_tag that was introduced upon this method execution.
        By default, the files on disk will be also renamed to reflect the
        genome configuration file changes

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of particular asset of interest
        :param str tag: name of the tag that identifies the asset of interest
        :param str new_tag: name of particular the new tag
        :param bool files: whether the asset files on disk should be renamed
        :raise ValueError: when the original tag is not specified
        :return bool: a logical indicating whether the tagging was successful
        """
        if any([c in new_tag for c in TAG_NAME_BANNED_CHARS]):
            raise ValueError(
                f"The tag name can't consist of characters: {TAG_NAME_BANNED_CHARS}"
            )
        self.run_plugins(PRE_TAG_HOOK)
        ori_path = self.seek_src(
            genome, asset, tag, enclosing_dir=True, strict_exists=True
        )
        alias_ori_path = self.seek(
            genome, asset, tag, enclosing_dir=True, strict_exists=True
        )
        new_path = os.path.abspath(os.path.join(ori_path, os.pardir, new_tag))
        if self.file_path:
            with self as r:
                if not r.cfg_tag_asset(genome, asset, tag, new_tag, force):
                    sys.exit(0)
        else:
            if not self.cfg_tag_asset(genome, asset, tag, new_tag, force):
                sys.exit(0)
        if not files:
            self.run_plugins(POST_TAG_HOOK)
            return
        try:
            if os.path.exists(new_path):
                _remove(new_path)
            os.rename(ori_path, new_path)
            _LOGGER.info("Renamed directory: {}".format(new_path))
            self._symlink_alias(genome, asset, new_tag)
            _remove(alias_ori_path)
        except FileNotFoundError:
            _LOGGER.warning(
                "Could not rename original asset tag directory '{}'"
                " to the new one '{}'".format(ori_path, new_path)
            )
        else:
            if self.file_path:
                with self as r:
                    r.cfg_remove_assets(genome, asset, tag, relationships=False)
            else:
                self.cfg_remove_assets(genome, asset, tag, relationships=False)
            _LOGGER.debug(
                "Asset '{}/{}' tagged with '{}' has been removed from"
                " the genome config".format(genome, asset, tag)
            )
            _LOGGER.debug(
                "Original asset has been moved from '{}' to '{}'".format(
                    ori_path, new_path
                )
            )
        self.run_plugins(POST_TAG_HOOK)

    def cfg_tag_asset(self, genome, asset, tag, new_tag, force=False):
        """
        Retags the asset selected by the tag with the new_tag.
        Prompts if default already exists and overrides upon confirmation.

        This method does not override the original asset entry in the
        RefGenConf object. It creates its copy and tags it with the new_tag.
        Additionally, if the retagged asset has any children their parent will
         be retagged as new_tag that was introduced upon this method execution.

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of particular asset of interest
        :param str tag: name of the tag that identifies the asset of interest
        :param str new_tag: name of particular the new tag
        :param bool force: force any actions that require approval
        :raise ValueError: when the original tag is not specified
        :return bool: a logical indicating whether the tagging was successful
        """
        self._assert_gat_exists(genome, asset, tag)
        asset_mapping = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset]
        if tag is None:
            ts = ", ".join(get_asset_tags(asset_mapping))
            raise ValueError(
                f"You must explicitly specify the tag of the asset"
                f" you want to reassign. Currently defined tags "
                f"for '{genome}/{asset}' are: {ts}"
            )
        if new_tag in asset_mapping[CFG_ASSET_TAGS_KEY]:
            if not force and not query_yes_no(
                f"You already have a '{asset}' asset tagged as "
                f"'{new_tag}', do you wish to override?"
            ):
                _LOGGER.info("Tag action aborted by the user")
                return
        children = []
        parents = []
        if CFG_ASSET_CHILDREN_KEY in asset_mapping[CFG_ASSET_TAGS_KEY][tag]:
            children = asset_mapping[CFG_ASSET_TAGS_KEY][tag][CFG_ASSET_CHILDREN_KEY]
        if CFG_ASSET_PARENTS_KEY in asset_mapping[CFG_ASSET_TAGS_KEY][tag]:
            parents = asset_mapping[CFG_ASSET_TAGS_KEY][tag][CFG_ASSET_PARENTS_KEY]
        if len(children) > 0 or len(parents) > 0:
            if not force and not query_yes_no(
                f"The asset '{genome}/{asset}:{tag}' has {len(children)} "
                f"children and {len(parents)} parents. Refgenie will update"
                f" the relationship data. Do you want to proceed?"
            ):
                _LOGGER.info("Tag action aborted by the user")
                return False
            # updates children's parents
            self._update_relatives_tags(
                genome, asset, tag, new_tag, children, update_children=False
            )
            # updates parents' children
            self._update_relatives_tags(
                genome, asset, tag, new_tag, parents, update_children=True
            )
        self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][CFG_ASSET_TAGS_KEY][
            new_tag
        ] = asset_mapping[CFG_ASSET_TAGS_KEY][tag]
        if (
            CFG_ASSET_DEFAULT_TAG_KEY in asset_mapping
            and asset_mapping[CFG_ASSET_DEFAULT_TAG_KEY] == tag
        ):
            self.set_default_pointer(
                genome, asset, new_tag, force_exists=True, force_fasta=True
            )
        self.cfg_remove_assets(genome, asset, tag)
        return True

    def _update_relatives_tags(
        self, genome, asset, tag, new_tag, relatives, update_children
    ):
        """
        Internal method used for tags updating in the 'asset_parents' section in the
        list of children.

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of particular asset of interest
        :param str tag: name of the tag that identifies the asset of interest
        :param str new_tag: name of particular the new tag
        :param list[str] relatives: relatives to be updated. Format: ["asset_name:tag",
            "asset_name1:tag1"]
        :param bool update_children: whether the children of the selected relatives
            should be updated.
        """
        relative_key = (
            CFG_ASSET_CHILDREN_KEY if update_children else CFG_ASSET_PARENTS_KEY
        )
        for r in relatives:
            _LOGGER.debug(
                "updating {} in '{}'".format(
                    "children" if update_children else "parents", r
                )
            )
            r_data = prp(r)
            try:
                self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][r_data["item"]][
                    CFG_ASSET_TAGS_KEY
                ][r_data["tag"]]
            except KeyError:
                _LOGGER.warning(
                    "The {} asset of '{}/{}' does not exist: {}".format(
                        "parent" if update_children else "child", genome, asset, r
                    )
                )
                continue
            updated_relatives = []
            if (
                relative_key
                in self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][r_data["item"]][
                    CFG_ASSET_TAGS_KEY
                ][r_data["tag"]]
            ):
                relatives = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][
                    r_data["item"]
                ][CFG_ASSET_TAGS_KEY][r_data["tag"]][relative_key]
                for relative in relatives:
                    ori_relative_data = prp(relative)
                    ori_relative_data["namespace"] = self.get_genome_alias_digest(
                        alias=ori_relative_data["namespace"], fallback=True
                    )
                    if (
                        ori_relative_data["item"] == asset
                        and ori_relative_data["tag"] == tag
                    ):
                        ori_relative_data["tag"] = new_tag
                        updated_relatives.append(
                            "{}/{}:{}".format(
                                ori_relative_data["namespace"], asset, new_tag
                            )
                        )
                    else:
                        updated_relatives.append(
                            "{}/{}:{}".format(
                                ori_relative_data["namespace"],
                                ori_relative_data["item"],
                                ori_relative_data["tag"],
                            )
                        )
            self.update_relatives_assets(
                genome,
                r_data["item"],
                r_data["tag"],
                updated_relatives,
                update_children,
            )
            self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][r_data["item"]][
                CFG_ASSET_TAGS_KEY
            ][r_data["tag"]][relative_key] = updated_relatives

    def pull(
        self,
        genome,
        asset,
        tag,
        unpack=True,
        force=None,
        force_large=None,
        size_cutoff=10,
        get_json_url=lambda server, operation_id: construct_request_url(
            server, operation_id
        ),
        build_signal_handler=_handle_sigint,
    ):
        """
        Download and possibly unpack one or more assets for a given ref gen.

        :param str genome: name of a reference genome assembly of interest
        :param str asset: name of particular asset to fetch
        :param str tag: name of particular tag to fetch
        :param bool unpack: whether to unpack a tarball
        :param bool | NoneType force: how to handle case in which asset path
            already exists; null for prompt (on a per-asset basis), False to
            effectively auto-reply No to the prompt to replace existing file,
            and True to auto-replay Yes for existing asset replacement.
        :param bool | NoneType force_large: how to handle case in large (> 5GB)
            asset is to be pulled; null for prompt (on a per-asset basis), False
            to effectively auto-reply No to the prompt,
            and True to auto-replay Yes
        :param float size_cutoff: maximum archive file size to download with
            no prompt
        :param function(str, str) -> str get_json_url: how to build URL from
            genome server URL base, genome, and asset
        :param function(str) -> function build_signal_handler: how to create
            a signal handler to use during the download; the single argument
            to this function factory is the download filepath
        :return (list[str], dict, str): a list of genome, asset, tag names
            and a key-value pair with which genome config file should be updated
            if pull succeeds, else asset key and a null value
        :raise refgenconf.UnboundEnvironmentVariablesError: if genome folder
            path contains any env. var. that's unbound
        :raise refgenconf.RefGenConfError: if the object update is requested in
            a non-writable state
        """
        self.run_plugins(PRE_PULL_HOOK)

        def _null_return():
            self.run_plugins(POST_PULL_HOOK)
            return gat, None, None

        def _raise_unpack_error():
            raise NotImplementedError(
                "Option to not extract tarballs is not yet supported."
            )

        num_servers = 0
        bad_servers = []
        no_asset_json = []
        alias = genome
        gat = [genome, asset, tag]
        if CFG_SERVERS_KEY not in self or self[CFG_SERVERS_KEY] is None:
            _LOGGER.error("You are not subscribed to any asset servers")
            return _null_return()

        good_servers = [
            s for s in self[CFG_SERVERS_KEY] if get_json_url(s, API_ID_DIGEST)
        ]

        _LOGGER.info(f"Compatible refgenieserver instances: {good_servers}")

        for server_url in good_servers:
            try:
                genome = self.get_genome_alias_digest(alias=alias)
            except yacman.UndefinedAliasError:
                _LOGGER.info(f"No local digest for genome alias: {genome}")
                if not self.set_genome_alias(
                    genome=alias, servers=[server_url], create_genome=True
                ):
                    continue
                genome = self.get_genome_alias_digest(alias=alias)

            num_servers += 1
            try:
                determined_tag = (
                    send_data_request(
                        get_json_url(server_url, API_ID_DEFAULT_TAG).format(
                            genome=genome, asset=asset
                        )
                    )
                    if tag is None
                    else tag
                )
            except DownloadJsonError as e:
                _LOGGER.warning(
                    f"Could not retrieve tag from: {server_url}. Caught exception: {e}"
                )
                bad_servers.append(server_url)
                continue
            else:
                determined_tag = str(determined_tag)
                _LOGGER.debug(f"Determined tag: {determined_tag}")
                unpack or _raise_unpack_error()
            gat = [genome, asset, determined_tag]
            url_asset_attrs = get_json_url(server_url, API_ID_ASSET_ATTRS).format(
                genome=genome, asset=asset
            )
            url_genome_attrs = get_json_url(server_url, API_ID_GENOME_ATTRS).format(
                genome=genome
            )
            url_archive = get_json_url(server_url, API_ID_ARCHIVE).format(
                genome=genome, asset=asset
            )

            try:
                archive_data = send_data_request(
                    url_asset_attrs, params={"tag": determined_tag}
                )
            except DownloadJsonError:
                no_asset_json.append(server_url)
                if num_servers == len(good_servers):
                    _LOGGER.error(
                        f"'{genome}/{asset}:{determined_tag}' not "
                        f"available on any of the following servers: "
                        f"{', '.join(self[CFG_SERVERS_KEY])}"
                    )
                    return _null_return()
                continue
            else:
                _LOGGER.debug("Determined server URL: {}".format(server_url))
                genome_archive_data = send_data_request(url_genome_attrs)

            if sys.version_info[0] == 2:
                archive_data = asciify_json_dict(archive_data)

            # local directory that the asset data will be stored in
            tag_dir = os.path.dirname(self.filepath(*gat))
            # local target path for the saved archive
            tardir = os.path.join(self.data_dir, genome, asset)
            tarpath = os.path.join(tardir, asset + "__" + determined_tag + ".tgz")
            # check if the genome/asset:tag exists and get request user decision
            if os.path.exists(tag_dir):

                def preserve():
                    _LOGGER.info(f"Preserving existing: {tag_dir}")
                    return _null_return()

                if force is False:
                    return preserve()
                elif force is None:
                    if not query_yes_no(f"Replace existing ({tag_dir})?", "no"):
                        return preserve()
                    else:
                        _LOGGER.debug(f"Overwriting: {tag_dir}")
                else:
                    _LOGGER.debug(f"Overwriting: {tag_dir}")

            # check asset digests local-server match for each parent
            [
                self._chk_digest_if_avail(
                    genome, x, archive_data[CFG_ASSET_CHECKSUM_KEY]
                )
                for x in archive_data[CFG_ASSET_PARENTS_KEY]
                if CFG_ASSET_PARENTS_KEY in archive_data
            ]

            bundle_name = "{}/{}:{}".format(*gat)
            archsize = archive_data[CFG_ARCHIVE_SIZE_KEY]
            _LOGGER.debug(f"'{bundle_name}' archive size: {archsize}")

            if not force_large and _is_large_archive(archsize, size_cutoff):
                if force_large is False:
                    _LOGGER.info(
                        "Skipping pull of {}/{}:{}; size: {}".format(*gat, archsize)
                    )
                    return _null_return()
                if not query_yes_no(
                    "This archive exceeds the size cutoff ({} > {:.1f}GB). "
                    "Do you want to proceed?".format(archsize, size_cutoff)
                ):
                    _LOGGER.info(
                        "Skipping pull of {}/{}:{}; size: {}".format(*gat, archsize)
                    )
                    return _null_return()

            if not os.path.exists(tardir):
                _LOGGER.debug(f"Creating directory: {tardir}")
                os.makedirs(tardir)

            # Download the file from `url` and save it locally under `filepath`:
            _LOGGER.info(f"Downloading URL: {url_archive}")
            try:
                signal.signal(signal.SIGINT, build_signal_handler(tarpath))
                _download_url_progress(
                    url_archive, tarpath, bundle_name, params={"tag": determined_tag}
                )
            except HTTPError:
                _LOGGER.error(
                    "Asset archive '{}/{}:{}' is missing on the "
                    "server: {s}".format(*gat, s=server_url)
                )
                if server_url == self[CFG_SERVERS_KEY][-1]:
                    # it this was the last server on the list, return
                    return _null_return()
                else:
                    _LOGGER.info("Trying next server")
                    # set the tag value back to what user requested
                    determined_tag = tag
                    continue
            except ConnectionRefusedError as e:
                _LOGGER.error(str(e))
                _LOGGER.error(
                    f"Server {server_url}/{API_VERSION} refused "
                    f"download. Check your internet settings"
                )
                return _null_return()
            except ContentTooShortError as e:
                _LOGGER.error(str(e))
                _LOGGER.error(f"'{bundle_name}' download incomplete")
                return _null_return()
            else:
                _LOGGER.info(f"Download complete: {tarpath}")

            new_checksum = checksum(tarpath)
            old_checksum = archive_data and archive_data.get(CFG_ARCHIVE_CHECKSUM_KEY)
            if old_checksum and new_checksum != old_checksum:
                _LOGGER.error(
                    f"Downloaded archive ('{tarpath}') checksum "
                    f"mismatch: ({new_checksum}, {old_checksum})"
                )
                return _null_return()
            else:
                _LOGGER.debug(f"Matched checksum: '{old_checksum}'")
            # successfully downloaded tarball; untar it
            if unpack and tarpath.endswith(".tgz"):
                _LOGGER.info(f"Extracting asset tarball: {tarpath}")
                untar(tarpath, tardir)
                os.remove(tarpath)

            if self.file_path:
                with self as rgc:
                    [
                        rgc.chk_digest_update_child(
                            gat[0], x, "{}/{}:{}".format(*gat), server_url
                        )
                        for x in archive_data[CFG_ASSET_PARENTS_KEY]
                        if CFG_ASSET_PARENTS_KEY in archive_data
                    ]
                    rgc.update_tags(
                        *gat,
                        data={
                            attr: archive_data[attr]
                            for attr in ATTRS_COPY_PULL
                            if attr in archive_data
                        },
                    )
                    rgc.set_default_pointer(*gat)
                    rgc.update_genomes(genome=genome, data=genome_archive_data)
            else:
                [
                    self.chk_digest_update_child(
                        gat[0], x, "{}/{}:{}".format(*gat), server_url
                    )
                    for x in archive_data[CFG_ASSET_PARENTS_KEY]
                    if CFG_ASSET_PARENTS_KEY in archive_data
                ]
                self.update_tags(
                    *gat,
                    data={
                        attr: archive_data[attr]
                        for attr in ATTRS_COPY_PULL
                        if attr in archive_data
                    },
                )
                self.set_default_pointer(*gat)
                self.update_genomes(genome=genome, data=genome_archive_data)
            if asset == "fasta":
                self.initialize_genome(
                    fasta_path=self.seek_src(*gat), alias=alias, fasta_unzipped=True
                )
            self.run_plugins(POST_PULL_HOOK)
            self._symlink_alias(*gat)
            return gat, archive_data, server_url

    def get_genome_alias_digest(self, alias, fallback=False):
        """
        Get the human readable alias for a genome digest

        :param str alias: alias to find digest for
        :param bool fallback: whether to return the query alias in case
            of failure and in case it is one of the digests
        :return str: genome digest
        :raise UndefinedAliasError: if the specified alias has been assigned to
            any digests
        """
        try:
            return self[CFG_GENOMES_KEY].get_key(alias=alias)
        except (yacman.UndefinedAliasError, AttributeError):
            if not fallback:
                raise
            if alias in self.genome_aliases.values():
                return alias
            raise

    def get_genome_alias(self, digest, fallback=False, all_aliases=False):
        """
        Get the human readable alias for a genome digest

        :param str digest: digest to find human-readable alias for
        :param bool fallback: whether to return the query digest in case
            of failure
        :param bool all_aliases: whether to return all aliases instead of just
            the first one
        :return str | list[str]: human-readable aliases
        :raise GenomeConfigFormatError: if "genome_digests" section does
            not exist in the config
        :raise UndefinedAliasError: if a no alias has been defined for the
            requested digest
        """
        try:
            res = self[CFG_GENOMES_KEY].get_aliases(key=digest)
            return res if all_aliases else res[0]
        except (yacman.UndefinedAliasError, AttributeError):
            if not fallback:
                raise
            if digest in self.genome_aliases.keys():
                return digest
            raise

    def remove_genome_aliases(self, digest, aliases=None):
        """
        Remove alias for a specified genome digest. This method will remove the
        digest both from the genomes object and from the aliases mapping
        in tbe config

        :param str digest: genome digest to remove an alias for
        :param list[str] aliases: a collection to aliases to remove for the
            genome. If not provided, all aliases for the digest will be remove
        :return bool: whether the removal has been performed
        """

        def _check_and_remove_alias(rgc, d, a):
            """
            Remove genome alias only if the alias can be remove successfully and
            genome exists
            """
            if rgc[CFG_GENOMES_KEY]:
                rmd = rgc[CFG_GENOMES_KEY].remove_aliases(key=d, aliases=a)
                if not rmd:
                    return rmd
                try:
                    rgc[CFG_GENOMES_KEY][d][CFG_ALIASES_KEY] = rgc[
                        CFG_GENOMES_KEY
                    ].get_aliases(d)
                except KeyError:
                    return []
                except yacman.UndefinedAliasError:
                    rgc[CFG_GENOMES_KEY][d][CFG_ALIASES_KEY] = []
                return rmd

        # get the symlink mapping before the removal for _remove_symlink_alias
        symlink_mapping = self.get_symlink_paths(genome=digest, all_aliases=True)
        if self.file_path:
            with self as r:
                removed_aliases = _check_and_remove_alias(r, digest, aliases)
        else:
            removed_aliases = _check_and_remove_alias(self, digest, aliases)
        if not removed_aliases:
            return [], []
        self._remove_symlink_alias(symlink_mapping, removed_aliases)
        return removed_aliases

    def set_genome_alias(
        self,
        genome,
        digest=None,
        servers=None,
        overwrite=False,
        reset_digest=False,
        create_genome=False,
        no_write=False,
        get_json_url=lambda server: construct_request_url(server, API_ID_ALIAS_DIGEST),
    ):
        """
        Assign a human-readable alias to a genome identifier.

        Genomes are identified by a unique identifier which is derived from the
        FASTA file (part of fasta asset). This way we can ensure genome
        provenance and compatibility with the server. This function maps a
        human-readable identifier to make referring to the genomes easier.

        :param str genome: name of the genome to assign to an identifier
        :param str digest: identifier to use
        :param bool overwrite: whether all the previously set aliases should be
            removed and just the current one stored
        :param bool no_write: whether to skip writing the alias to the file
        :return bool: whether the alias has been established
        """

        def _check_and_set_alias(rgc, d, a, create=False):
            """
            Set genome alias only if the key alias can be set successfully and
            genome exists or genome creation is forced
            """
            try:
                _assert_gat_exists(rgc[CFG_GENOMES_KEY], gname=digest)
            except MissingGenomeError:
                if not create:
                    raise
                rgc[CFG_GENOMES_KEY][d] = PXAM()

            sa, ra = rgc[CFG_GENOMES_KEY].set_aliases(
                aliases=a, key=d, overwrite=overwrite, reset_key=reset_digest
            )
            try:
                rgc[CFG_GENOMES_KEY][d][CFG_ALIASES_KEY] = rgc[
                    CFG_GENOMES_KEY
                ].get_aliases(d)
            except KeyError:
                return [], []
            _LOGGER.info(
                f"Set genome alias ({d}: {', '.join(a) if isinstance(a, list) else a})"
            )
            return sa, ra

        if not digest:
            if isinstance(genome, list):
                if len(genome) > 1:
                    raise NotImplementedError("Can look up just one digest at a time")
                else:
                    genome = genome[0]
            cnt = 0
            if servers is None:
                servers = self[CFG_SERVERS_KEY]
            for server in servers:
                cnt += 1
                url_alias_template = get_json_url(server=server)
                if url_alias_template is None:
                    continue
                url_alias = url_alias_template.format(alias=genome)
                _LOGGER.info(f"Setting '{genome}' identity with server: {url_alias}")
                try:
                    digest = send_data_request(url_alias)
                except DownloadJsonError:
                    if cnt == len(servers):
                        _LOGGER.error(
                            f"Genome '{genome}' not available on any of the "
                            f"following servers: {', '.join(servers)}"
                        )
                        return False
                    continue
                _LOGGER.info(f"Determined digest for local '{genome}' alias: {digest}")
                break

        # get the symlink mapping before the removal for _remove_symlink_alias
        symlink_mapping = self.get_symlink_paths(genome=digest, all_aliases=True)
        if self.file_path and not no_write:
            with self as r:
                set_aliases, removed_aliases = _check_and_set_alias(
                    rgc=r, d=digest, a=genome, create=create_genome
                )
            self._remove_symlink_alias(symlink_mapping, removed_aliases)
            self._symlink_alias(genome=digest)
        else:
            set_aliases, removed_aliases = _check_and_set_alias(
                rgc=self, d=digest, a=genome, create=create_genome
            )
        if not set_aliases:
            return False
        return True

    def initialize_genome(
        self, fasta_path, alias, fasta_unzipped=False, skip_alias_write=False
    ):
        """
        Initialize a genome

        Create a JSON file with Annotated Sequence Digests (ASDs)
        for the FASTA file in the genome directory.

        :param str fasta_path: path to a FASTA file to initialize genome with
        :param str alias: alias to set for the genome
        :param bool skip_alias_write: whether to skip writing the alias to the file
        :return str, list[dict[]]: human-readable name for the genome
        """
        _LOGGER.info("Initializing genome: {}".format(alias))
        if not os.path.isfile(fasta_path):
            raise FileNotFoundError(
                "Can't initialize genome; FASTA file does "
                "not exist: {}".format(fasta_path)
            )
        ssc = SeqColClient({})
        d, _ = ssc.load_fasta(fasta_path, gzipped=not fasta_unzipped)
        # retrieve annotated sequence digests list to save in a JSON file
        asdl = ssc.retrieve(druid=d)
        pth = self.get_asds_path(d)
        os.makedirs(os.path.dirname(pth), exist_ok=True)
        with open(pth, "w") as jfp:
            json.dump(asdl, jfp)
        _LOGGER.debug("Saved ASDs to JSON: {}".format(pth))
        self.set_genome_alias(
            genome=alias,
            digest=d,
            overwrite=True,
            create_genome=True,
            no_write=skip_alias_write,
        )
        return d, asdl

    def get_asds_path(self, genome):
        """
        Get path to the Annotated Sequence Digests JSON file for a given genome.
        Note that the path and/or genome may not exist.

        :param str genome: genome name
        :return str: ASDs path
        """
        return os.path.join(self.data_dir, genome, f"{genome}__ASDs.json")

    def remove_asset_from_relatives(self, genome, asset, tag):
        """
        Remove any relationship links associated with the selected asset

        :param str genome: genome to be removed from its relatives' relatives list
        :param str asset: asset to be removed from its relatives' relatives list
        :param str tag: tag to be removed from its relatives' relatives list
        """
        to_remove = "{}/{}:{}".format(
            self.get_genome_alias_digest(alias=genome, fallback=True), asset, tag
        )
        for rel_type in CFG_ASSET_RELATIVES_KEYS:
            tmp = CFG_ASSET_RELATIVES_KEYS[
                len(CFG_ASSET_RELATIVES_KEYS)
                - 1
                - CFG_ASSET_RELATIVES_KEYS.index(rel_type)
            ]
            tag_data = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                CFG_ASSET_TAGS_KEY
            ][tag]
            if rel_type not in tag_data:
                continue
            for rel in tag_data[rel_type]:
                parsed = prp(rel)
                _LOGGER.debug("Removing '{}' from '{}' {}".format(to_remove, rel, tmp))
                try:
                    self[CFG_GENOMES_KEY][parsed["namespace"] or genome][
                        CFG_ASSETS_KEY
                    ][parsed["item"]][CFG_ASSET_TAGS_KEY][parsed["tag"]][tmp].remove(
                        to_remove
                    )
                except (KeyError, ValueError):
                    pass

    def update_relatives_assets(
        self, genome, asset, tag=None, data=None, children=False
    ):
        """
        A convenience method which wraps the update assets and uses it to update the
        asset relatives of an asset.

        :param str genome: genome to be added/updated
        :param str asset: asset to be added/updated
        :param str tag: tag to be added/updated
        :param list data: asset parents or children to be added/updated
        :param bool children: a logical indicating whether the relationship to be
            added is 'children'
        :return RefGenConf: updated object
        """
        tag = tag or self.get_default_tag(genome, asset)
        relationship = CFG_ASSET_CHILDREN_KEY if children else CFG_ASSET_PARENTS_KEY
        if _check_insert_data(data, list, "data"):
            # creates/asserts the genome/asset:tag combination
            self.update_tags(genome, asset, tag)
            tag_data = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                CFG_ASSET_TAGS_KEY
            ][tag]
            tag_data.setdefault(relationship, list())
            tag_data[relationship] = _extend_unique(
                tag_data[relationship],
                data,
            )

    def update_seek_keys(self, genome, asset, tag=None, keys=None, force_digest=None):
        """
        A convenience method which wraps the updated assets and uses it to
        update the seek keys for a tagged asset.

        :param str genome: genome to be added/updated
        :param str asset: asset to be added/updated
        :param str tag: tag to be added/updated
        :param str force_digest: digest to force update of. The alias will
            not be converted to the digest, even if provided.
        :param Mapping keys: seek_keys to be added/updated
        :return RefGenConf: updated object
        """
        tag = tag or self.get_default_tag(genome, asset)
        if _check_insert_data(keys, Mapping, "keys"):
            self.update_tags(genome, asset, tag, force_digest=force_digest)
            asset = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset]
            _safe_setdef(asset[CFG_ASSET_TAGS_KEY][tag], CFG_SEEK_KEYS_KEY, PXAM())
            asset[CFG_ASSET_TAGS_KEY][tag][CFG_SEEK_KEYS_KEY].update(keys)
        return self

    def update_tags(self, genome, asset=None, tag=None, data=None, force_digest=None):
        """
        Updates the genomes in RefGenConf object at any level.
        If a requested genome-asset-tag mapping is missing, it will be created

        :param str genome: genome to be added/updated
        :param str asset: asset to be added/updated
        :param str tag: tag to be added/updated
        :param str force_digest: digest to force update of. The alias will
            not be converted to the digest, even if provided.
        :param Mapping data: data to be added/updated
        :return RefGenConf: updated object
        """
        if _check_insert_data(genome, str, "genome"):
            genome = force_digest or self.get_genome_alias_digest(
                alias=genome, fallback=True
            )
            _safe_setdef(self[CFG_GENOMES_KEY], genome, PXAM())
            if _check_insert_data(asset, str, "asset"):
                _safe_setdef(self[CFG_GENOMES_KEY][genome], CFG_ASSETS_KEY, PXAM())
                _safe_setdef(
                    self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY], asset, PXAM()
                )
                if _check_insert_data(tag, str, "tag"):
                    _safe_setdef(
                        self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset],
                        CFG_ASSET_TAGS_KEY,
                        PXAM(),
                    )
                    _safe_setdef(
                        self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                            CFG_ASSET_TAGS_KEY
                        ],
                        tag,
                        PXAM(),
                    )
                    if _check_insert_data(data, Mapping, "data"):
                        self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                            CFG_ASSET_TAGS_KEY
                        ][tag].update(data)
        return self

    def update_assets(self, genome, asset=None, data=None, force_digest=None):
        """
        Updates the genomes in RefGenConf object at any level.
        If a requested genome-asset mapping is missing, it will be created

        :param str genome: genome to be added/updated
        :param str asset: asset to be added/updated
        :param str force_digest: digest to force update of. The alias will
            not be converted to the digest, even if provided.
        :param Mapping data: data to be added/updated
        :return RefGenConf: updated object
        """
        if _check_insert_data(genome, str, "genome"):
            genome = force_digest or self.get_genome_alias_digest(
                alias=genome, fallback=True
            )
            _safe_setdef(self[CFG_GENOMES_KEY], genome, PXAM())
            if _check_insert_data(asset, str, "asset"):
                _safe_setdef(self[CFG_GENOMES_KEY][genome], CFG_ASSETS_KEY, PXAM())
                _safe_setdef(
                    self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY], asset, PXAM()
                )
                if _check_insert_data(data, Mapping, "data"):
                    self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset].update(data)
        return self

    def remove(
        self, genome, asset, tag=None, relationships=True, files=True, force=False
    ):
        """
        Remove data associated with a specified genome:asset:tag combination.
        If no tags are specified, the entire asset is removed from the genome.

        If no more tags are defined for the selected genome:asset after tag removal,
        the parent asset will be removed as well
        If no more assets are defined for the selected genome after asset removal,
        the parent genome will be removed as well

        :param str genome: genome to be removed
        :param str asset: asset package to be removed
        :param str tag: tag to be removed
        :param bool relationships: whether the asset being removed should
            be removed from its relatives as well
        :param bool files: whether the asset files from disk should be removed
        :param bool force: whether the removal prompts should be skipped
        :raise TypeError: if genome argument type is not a list or str
        :return RefGenConf: updated object
        """
        tag = tag or self.get_default_tag(genome, asset, use_existing=False)
        if files:
            req_dict = {
                "genome": self.get_genome_alias_digest(genome, fallback=True),
                "asset": asset,
                "tag": tag,
            }
            _LOGGER.debug("Attempting removal: {}".format(req_dict))
            if not force and not query_yes_no(
                "Remove '{}/{}:{}'?".format(genome, asset, tag)
            ):
                _LOGGER.info("Action aborted by the user")
                return
            removed = []
            asset_path = self.seek_src(
                genome, asset, tag, enclosing_dir=True, strict_exists=False
            )
            alias_asset_paths = self.seek(
                genome,
                asset,
                tag,
                enclosing_dir=True,
                strict_exists=False,
                all_aliases=True,
            )
            if os.path.exists(asset_path):
                removed.append(_remove(asset_path))
                removed.extend([_remove(p) for p in alias_asset_paths])
                if self.file_path:
                    with self as r:
                        r.cfg_remove_assets(genome, asset, tag, relationships)
                else:
                    self.cfg_remove_assets(genome, asset, tag, relationships)
            else:
                _LOGGER.warning(
                    "Selected asset does not exist on disk ({}). "
                    "Removing from genome config.".format(asset_path)
                )
                if self.file_path:
                    with self as r:
                        r.cfg_remove_assets(genome, asset, tag, relationships)
                        return
                else:
                    self.cfg_remove_assets(genome, asset, tag, relationships)
                    return
            try:
                self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset]
            except (KeyError, TypeError):
                asset_dir = os.path.abspath(os.path.join(asset_path, os.path.pardir))
                alias_asset_dirs = [
                    os.path.abspath(os.path.join(p, os.path.pardir))
                    for p in alias_asset_paths
                ]
                _entity_dir_removal_log(asset_dir, "asset", req_dict, removed)
                removed.extend([_remove(p) for p in alias_asset_dirs])
                try:
                    self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY]
                except (KeyError, TypeError):
                    genome_dir = os.path.abspath(
                        os.path.join(asset_dir, os.path.pardir)
                    )
                    alias_genome_dirs = [
                        os.path.abspath(os.path.join(p, os.path.pardir))
                        for p in alias_asset_dirs
                    ]
                    _entity_dir_removal_log(genome_dir, "genome", req_dict, removed)
                    removed.extend([_remove(p) for p in alias_genome_dirs])
                    try:
                        if self.file_path:
                            with self as r:
                                del r[CFG_GENOMES_KEY][genome]
                        else:
                            del self[CFG_GENOMES_KEY][genome]
                    except (KeyError, TypeError):
                        _LOGGER.debug(
                            "Could not remove genome '{}' from the config; it "
                            "does not exist".format(genome)
                        )
            _LOGGER.info(f"Successfully removed entities:{block_iter_repr(removed)}")
        else:
            if self.file_path:
                with self as r:
                    r.cfg_remove_assets(genome, asset, tag, relationships)
            else:
                self.cfg_remove_assets(genome, asset, tag, relationships)

    def cfg_remove_assets(self, genome, asset, tag=None, relationships=True):
        """
        Remove data associated with a specified genome:asset:tag combination.
        If no tags are specified, the entire asset is removed from the genome.

        If no more tags are defined for the selected genome:asset after tag removal,
        the parent asset will be removed as well
        If no more assets are defined for the selected genome after asset removal,
        the parent genome will be removed as well

        :param str genome: genome to be removed
        :param str asset: asset package to be removed
        :param str tag: tag to be removed
        :param bool relationships: whether the asset being removed should
            be removed from its relatives as well
        :raise TypeError: if genome argument type is not a list or str
        :return RefGenConf: updated object
        """

        def _del_if_empty(obj, attr, alt=None):
            """
            Internal function for Mapping attribute deleting.
            Check if attribute exists and delete it if its length is zero.

            :param Mapping obj: an object to check
            :param str attr: Mapping attribute of interest
            :param list[Mapping, str] alt: a list of length 2 that indicates alternative
            Mapping-attribute combination to remove
            """
            if attr in obj and len(obj[attr]) == 0:
                if alt is None:
                    del obj[attr]
                else:
                    if alt[1] in alt[0]:
                        del alt[0][alt[1]]

        tag = tag or self.get_default_tag(genome, asset)
        if _check_insert_data(genome, str, "genome"):
            if _check_insert_data(asset, str, "asset"):
                if _check_insert_data(tag, str, "tag"):
                    if relationships:
                        self.remove_asset_from_relatives(genome, asset, tag)
                    del self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                        CFG_ASSET_TAGS_KEY
                    ][tag]
                    _del_if_empty(
                        self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset],
                        CFG_ASSET_TAGS_KEY,
                        [self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY], asset],
                    )
                    _del_if_empty(self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY], asset)
                    _del_if_empty(
                        self[CFG_GENOMES_KEY][genome],
                        CFG_ASSETS_KEY,
                        [self[CFG_GENOMES_KEY], genome],
                    )
                    _del_if_empty(self[CFG_GENOMES_KEY], genome)
                    try:
                        default_tag = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][
                            asset
                        ][CFG_ASSET_DEFAULT_TAG_KEY]
                    except KeyError:
                        pass
                    else:
                        if default_tag == tag:
                            del self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                                CFG_ASSET_DEFAULT_TAG_KEY
                            ]
                    if len(self[CFG_GENOMES_KEY]) == 0:
                        self[CFG_GENOMES_KEY] = None
        return self

    def update_genomes(self, genome, data=None, force_digest=None):
        """
        Updates the genomes in RefGenConf object at any level.
        If a requested genome is missing, it will be added

        :param str genome: genome to be added/updated
        :param str force_digest: digest to force update of. The alias will
            not be converted to the digest, even if provided.
        :param Mapping data: data to be added/updated
        :return RefGenConf: updated object
        """
        if _check_insert_data(genome, str, "genome"):
            genome = force_digest or self.get_genome_alias_digest(
                alias=genome, fallback=True
            )
            _safe_setdef(self[CFG_GENOMES_KEY], genome, PXAM({CFG_ASSETS_KEY: PXAM()}))
            if _check_insert_data(data, Mapping, "data"):
                self[CFG_GENOMES_KEY][genome].update(data)
        return self

    def _update_genome_servers(self, url, reset=False):
        """
        Update the list of genome_servers.

        Use reset argument to overwrite the current list. Otherwise the current
        one will be appended to.

        :param list[str] | str url: url(s) to update the genome_servers list with
        :param bool reset: whether the current list should be overwritten
        """
        if CFG_SERVERS_KEY in self:
            self[CFG_SERVERS_KEY] = _extend_unique(
                [] if reset else self[CFG_SERVERS_KEY], _make_list_of_str(url)
            )
        else:
            raise GenomeConfigFormatError(
                "The '{}' is missing. Can't update the server list".format(
                    CFG_SERVERS_KEY
                )
            )

    def subscribe(self, urls, reset=False, no_write=False):
        """
        Add URLs the list of genome_servers.

        Use reset argument to overwrite the current list.
        Otherwise the current one will be appended to.

        :param list[str] | str urls: urls to update the genome_servers list with
        :param bool reset: whether the current list should be overwritten
        """
        if self.file_path and not no_write:
            with self as r:
                r._update_genome_servers(url=urls, reset=reset)
        else:
            self._update_genome_servers(url=urls, reset=reset)
        _LOGGER.info(f"Subscribed to: {', '.join(urls)}")

    def unsubscribe(self, urls, no_write=False):
        """
        Remove URLs the list of genome_servers.

        :param list[str] | str urls: urls to update the genome_servers list with
        """
        unsub_list = []
        ori_servers = self[CFG_SERVERS_KEY]
        for s in urls:
            try:
                ori_servers.remove(s)
                unsub_list.append(s)
            except ValueError:
                _LOGGER.warning(
                    "URL '{}' not in genome_servers list: {}".format(s, ori_servers)
                )
        if self.file_path and not no_write:
            with self as r:
                r._update_genome_servers(ori_servers, reset=True)
        else:
            self._update_genome_servers(ori_servers, reset=True)
        if unsub_list:
            _LOGGER.info("Unsubscribed from: {}".format(", ".join(unsub_list)))

    def getseq(self, genome, locus, as_str=False):
        """
        Return the sequence found in a selected range and chromosome.
        Something like the refget protocol.

        :param str genome: name of the sequence identifier
        :param str locus: coordinates of desired sequence, e.g. 'chr1:1-10'
        :param bool as_str: whether to convert the resurned object to string
            and return just the sequence
        :return str | pyfaidx.FastaRecord | pyfaidx.Sequence: selected sequence
        """
        import pyfaidx

        fa = pyfaidx.Fasta(self.seek_src(genome, "fasta", strict_exists=True))
        locus_split = locus.split(":")
        chr = fa[locus_split[0]]
        if len(locus_split) == 1:
            return str(chr) if as_str else chr
        start, end = locus_split[1].split("-")
        _LOGGER.debug(
            "chr: '{}', start: '{}', end: '{}'".format(locus_split[0], start, end)
        )
        return str(chr[int(start) : int(end)]) if as_str else chr[int(start) : int(end)]

    def get_genome_attributes(self, genome):
        """
        Get the dictionary attributes, like checksum, contents, description.
        Does not return the assets.

        :param str genome: genome to get the attributes dict for
        :return Mapping[str, str]: available genome attributes
        """
        return {
            k: self[CFG_GENOMES_KEY][genome][k]
            for k in CFG_GENOME_ATTRS_KEYS
            if k in self[CFG_GENOMES_KEY][genome]
        }

    def is_asset_complete(self, genome, asset, tag):
        """
        Check whether all required tag attributes are defined in the RefGenConf object.
        This is the way we determine tag completeness.

        :param str genome: genome to be checked
        :param str asset: asset package to be checked
        :param str tag: tag to be checked
        :return bool: the decision
        """
        tag_data = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
            CFG_ASSET_TAGS_KEY
        ][tag]
        return all([r in tag_data for r in REQ_TAG_ATTRS])

    def _invert_genomes(self, order=None):
        """Map each asset type/kind/name to a collection of assemblies.

        A configuration file encodes assets by genome, but in some use cases
        it's helpful to invert the direction of this mapping. The value of the
        asset key/name may differ by genome, so that information is
        necessarily lost in this inversion, but we can collect genome IDs by
        asset ID.

        :param function(str) -> object order: how to key genome IDs and asset
            names for sort
        :return OrderedDict[str, Iterable[str]] binding between asset kind/key/name
            and collection of reference genome assembly names for which the
            asset type is available
        """
        genomes = {}
        for g, am in self[CFG_GENOMES_KEY].items():
            for a in am[CFG_ASSETS_KEY].keys():
                genomes.setdefault(a, []).append(g)
        assets = sorted(genomes.keys(), key=order)
        return OrderedDict([(a, sorted(genomes[a], key=order)) for a in assets])

    def _chk_digest_if_avail(self, genome, remote_asset_name, remote_digest):
        """
        Check local asset digest against the remote one and populate children of the
        asset with the provided asset:tag.

        In case the local asset does not exist, the config is populated with the
        remote asset digest and children data

        :param str genome: name of the genome to check the asset digests for
        :param str remote_asset_name: asset and tag names, formatted like: asset:tag
        :param str server_url: address of the server to query for the digests
        :raise RefgenconfError: if the local digest does not match its remote counterpart
        """
        remote_asset_data = prp(remote_asset_name)
        asset = remote_asset_data["item"]
        tag = remote_asset_data["tag"]
        try:
            local_digest = self.id(genome, asset, tag)
            if remote_digest != local_digest:
                raise RemoteDigestMismatchError(asset, local_digest, remote_digest)
        except RefgenconfError:
            _LOGGER.debug(
                f"Could not find '{genome}/{asset}:{tag}' digest. Digest for this "
                f"parent will be populated with the server one after the pull"
            )
            return

    def chk_digest_update_child(
        self, genome, remote_asset_name, child_name, server_url
    ):
        """
        Check local asset digest against the remote one and populate children of the
        asset with the provided asset:tag.

        In case the local asset does not exist, the config is populated with the remote
         asset digest and children data

        :param str genome: name of the genome to check the asset digests for
        :param str remote_asset_name: asset and tag names, formatted like: asset:tag
        :param str child_name: name to be appended to the children of the parent
        :param str server_url: address of the server to query for the digests
        :raise RefgenconfError: if the local digest does not match its remote counterpart
        """
        remote_asset_data = prp(remote_asset_name)
        asset = remote_asset_data["item"]
        tag = remote_asset_data["tag"]
        asset_digest_url = construct_request_url(server_url, API_ID_DIGEST).format(
            genome=genome, asset=asset, tag=tag
        )
        try:
            remote_digest = send_data_request(asset_digest_url)
        except DownloadJsonError:
            return
        try:
            # we need to allow for missing seek_keys section so that the digest is
            # respected even from the previously populated 'incomplete asset' from
            # the server
            self._assert_gat_exists(
                genome,
                asset,
                tag,
                allow_incomplete=not self.is_asset_complete(genome, asset, tag),
            )
        except (KeyError, MissingAssetError, MissingGenomeError, MissingSeekKeyError):
            self.update_tags(
                genome, asset, tag, {CFG_ASSET_CHECKSUM_KEY: remote_digest}
            )
            _LOGGER.info(
                f"Could not find '{genome}/{asset}:{tag}' digest. "
                f"Populating with server data"
            )
        else:
            local_digest = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset][
                CFG_ASSET_TAGS_KEY
            ][tag][CFG_ASSET_CHECKSUM_KEY]
            if remote_digest != local_digest:
                raise RemoteDigestMismatchError(asset, local_digest, remote_digest)
        finally:
            self.update_relatives_assets(
                genome, asset, tag, [child_name], children=True
            )

    def id(self, genome, asset, tag=None):
        """
        Returns the digest for the specified asset.
        The defined default tag will be used if not provided as an argument

        :param str genome: genome identifier
        :param str asset: asset identifier
        :param str tag: tag identifier
        :return str: asset digest for the tag
        """
        self._assert_gat_exists(genome, asset, tag)
        tag = tag or self.get_default_tag(genome, asset)
        a = self[CFG_GENOMES_KEY][genome][CFG_ASSETS_KEY][asset]
        if CFG_ASSET_CHECKSUM_KEY in a[CFG_ASSET_TAGS_KEY][tag]:
            return a[CFG_ASSET_TAGS_KEY][tag][CFG_ASSET_CHECKSUM_KEY]
        raise MissingConfigDataError(
            "Digest does not exist for: {}/{}:{}".format(genome, asset, tag)
        )

    def compare(self, genome1, genome2, explain=False):
        """
        Check genomes compatibility level.
        Compares Annotated Sequence Digests (ASDs) -- digested sequences and
        metadata
        :param str genome1: name of the first genome to compare
        :param str genome2: name of the first genome to compare
        :param bool explain: whether the returned code explanation should
            be displayed
        :return int: compatibility code
        """

        def _get_asds_for_genome(rgc, genome):
            """
            Read JSON file containing ASDs for a specified genome
            :param refgenconf.RefGenConf rgc: object to find the genome for
            :param str genome: genome to find the file for
            :return list[dict]: list of ASDs, ready to compare
            """
            g = rgc.get_genome_alias(genome, fallback=True)
            error_msg = (
                f"File containing Annotated Sequence Digests (ASDs) not "
                f"found for genome: {g}. Must pull or build '{g}/fasta' again to "
                f"check the compatibility."
            )
            try:
                rgc.seek_src(genome, "fasta", strict_exists=True)
            except MissingSeekKeyError:
                raise MissingSeekKeyError(error_msg)
            json_file = rgc.get_asds_path(genome)
            if not os.path.exists(json_file):
                raise OSError(error_msg)
            with open(json_file, "r") as jfp:
                return json.load(jfp)

        return SeqColClient({}).compare_asds(
            _get_asds_for_genome(self, self.get_genome_alias_digest(genome1, True)),
            _get_asds_for_genome(self, self.get_genome_alias_digest(genome2, True)),
            explain=explain,
        )

    def populate(self, glob):
        """
        Populates *local* refgenie references from
        refgenie://genome/asset.seek_key:tag registry paths

        :param dict | str | list glob: String which may contain refgenie registry paths as
            values; or a dict, for which values may contain refgenie registry
            paths. Dict include nested dicts.
        :return dict | str | list: modified input dict with refgenie paths populated
        """
        return _populate_refgenie_registry_path(
            self, glob=glob, seek_method_name="seek"
        )

    def populater(self, glob, remote_class=None):
        """
        Populates *remote* refgenie references from
        refgenie://genome/asset:tag registry paths

        :param dict | str | list glob: String which may contain refgenie registry paths as
            values; or a dict, for which values may contain refgenie registry
            paths. Dict include nested dicts.
        :param str remote_class: remote data provider class, e.g. 'http' or 's3'
        :return dict | str | list: modified input dict with refgenie paths populated
        """
        return _populate_refgenie_registry_path(
            self,
            glob=glob,
            seek_method_name="seekr",
            remote_class=remote_class or "http",
        )

    def run_plugins(self, hook):
        """
        Runs all installed plugins for the specified hook.

        :param str hook: hook identifier
        """
        for name, func in self.plugins[hook].items():
            _LOGGER.debug("Running {} plugin: {}".format(hook, name))
            func(self)

    def write(self, filepath=None):
        """
        Write the contents to a file.
        If pre- and post-update plugins are defined, they will be executed automatically

        :param str filepath: a file path to write to
        :raise OSError: when the object has been created in a read only mode or other
            process has locked the file
        :raise TypeError: when the filepath cannot be determined.
            This takes place only if YacAttMap initialized with a Mapping as an input,
             not read from file.
        :raise OSError: when the write is called on an object with no write capabilities
            or when writing to a file that is locked by a different object
        :return str: the path to the created files
        """
        self.run_plugins(PRE_UPDATE_HOOK)
        try:
            path = super(RefGenConf, self).write(filepath=filepath, exclude_case=True)
        except ValidationError:
            _LOGGER.error("The changes were not written to the file")
            raise
        self.run_plugins(POST_UPDATE_HOOK)
        return path

    def _genome_asset_path(self, gname, aname, tname, seek_key, enclosing_dir):
        """
        Retrieve the raw path value for a particular asset for a particular genome.

        :param Mapping[str, Mapping[str, Mapping[str, object]]] genomes: nested
            collection of key-value pairs, keyed at top level on genome ID, then by
            asset name, then by asset attribute
        :param str gname: top level key to query -- genome ID, e.g. mm10
        :param str aname: second-level key to query -- asset name, e.g. fasta
        :param str tname: third-level key to query -- tag name, e.g. default
        :param str seek_key: fourth-level key to query -- tag name, e.g. chrom_sizes
        :param bool enclosing_dir: whether a path to the entire enclosing directory
            should be returned, e.g. for a fasta asset that has 3 seek_keys pointing
            to 3 files in an asset dir, that asset dir is returned
        :return str: raw path value for a particular asset for a particular genome
        :raise MissingGenomeError: if the given key-value pair collection does not
            contain as a top-level key the given genome ID
        :raise MissingAssetError: if the given key-value pair colelction does
            contain the given genome ID, but that key's mapping doesn't contain
            the given asset name as a key
        :raise GenomeConfigFormatError: if it's discovered during the query that
            the structure of the given genomes mapping suggests that it was
            parsed from an improperly formatted/structured genome config file.
        """
        self._assert_gat_exists(gname, aname, tname)
        asset_tag_data = self[CFG_GENOMES_KEY][gname][CFG_ASSETS_KEY][aname][
            CFG_ASSET_TAGS_KEY
        ][tname]
        if enclosing_dir:
            return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname)
        if seek_key is None:
            if aname in asset_tag_data[CFG_SEEK_KEYS_KEY]:
                seek_key = aname
            else:
                return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname)
        try:
            seek_key_value = asset_tag_data[CFG_SEEK_KEYS_KEY][seek_key]
            appendix = "" if seek_key_value == "." else seek_key_value
            return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname, appendix)
        except KeyError:
            raise MissingSeekKeyError(
                f"genome/asset:tag bundle '{gname}/{aname}:{tname}' exists, but "
                f"seek_key '{seek_key}' is missing"
            )

    def _get_genome_id(self, gname):
        """
        Get the actual genome name used in the object regardless of whether
        the query name is an actual name of an alias

        :param str gname: genome query name, can be the actual key or its alias
        :return str: genome id
        """
        self._assert_gat_exists(gname)
        if gname in self[CFG_GENOMES_KEY].keys():
            return gname
        return self[CFG_GENOMES_KEY].get_key(alias=gname)

    def _assert_gat_exists(self, gname, aname=None, tname=None, allow_incomplete=False):
        """
        Make sure the genome/asset:tag combination exists in the provided mapping and
        has any seek keys defined.
        Seek keys are required for the asset completeness.

        :param Mapping[str, Mapping[str, Mapping[str, object]]] genomes: nested
            collection of key-value pairs, keyed at top level on genome ID, then by
            asset name, then by asset attribute
        :param str gname: top level key to query -- genome ID, e.g. mm10
        :param str aname: second-level key to query -- asset name, e.g. fasta
        :param str tname: third-level key to query -- tag name, e.g. default
        :raise MissingGenomeError: if the given key-value pair collection does not
            contain as a top-level key the given genome ID
        :raise MissingAssetError: if the given key-value pair collection does
            contain the given genome ID, but that key's mapping doesn't contain
            the given asset name as a key
        :raise GenomeConfigFormatError: if it's discovered during the query that
            the structure of the given genomes mapping suggests that it was
            parsed from an improperly formatted/structured genome config file.
        """
        _LOGGER.debug(f"checking existence of: {gname}/{aname}:{tname}")
        try:
            genome = self[CFG_GENOMES_KEY][gname]
        except KeyError:
            raise MissingGenomeError(f"Your genomes do not include '{gname}'")
        if aname is not None:
            try:
                asset_data = genome[CFG_ASSETS_KEY][aname]
            except KeyError:
                raise MissingAssetError(
                    f"Genome '{gname}' exists, but asset '{aname}' is missing"
                )
            except TypeError:
                _raise_not_mapping(asset_data, "Asset section ")
            if tname is not None:
                try:
                    tag_data = asset_data[CFG_ASSET_TAGS_KEY][tname]
                except KeyError:
                    raise MissingTagError(
                        f"genome/asset bundle '{gname}/{aname}' exists, but tag "
                        f"'{tname}' is missing"
                    )
                except TypeError:
                    _raise_not_mapping(asset_data, "Asset section ")
                try:
                    tag_data[CFG_SEEK_KEYS_KEY]
                except KeyError:
                    if not allow_incomplete:
                        raise MissingSeekKeyError(
                            f"Asset incomplete. No seek keys are defined for "
                            f"'{gname}/{aname}:{tname}'. Build or pull the asset again."
                        )

    def _list_remote(
        self,
        url,
        genome,
    ):
        """
        List genomes and assets available remotely.

        :param url: location or ref genome config data
        :return str, str: text reps of remotely available genomes and assets
        """
        genomes_data = send_data_request(url, params={"includeSeekKeys": True})
        return (
            {g: data for g, data in genomes_data.items() if g in genome}
            if genome is not None
            else genomes_data
        )

    def _select_genomes(
        self, genome=None, strict=False, order=None, external_genomes=None
    ):
        """
        Safely select a subset of genomes

        :param list[str] | str genome: genomes that the assets should be found for
        :param bool strict: whether a non-existent genome should lead to a warning.
            Specific genome request is disregarded otherwise
        :param function(str) -> object order: a way to order the genomes in the output
        :param list external_genomes: a collection of genomes to use instead
            of the one defined in the object
        :raise TypeError: if genome argument type is not a list or str
        :return list: selected subset of genomes
        """
        if external_genomes:
            # expects remote genomes to be supplied as aliases; no conversion
            genomes = sorted(external_genomes, key=order)
        else:
            genomes = [
                self.get_genome_alias(x, fallback=True)
                for x in sorted(self[CFG_GENOMES_KEY].keys(), key=order)
            ]
        if not genome:
            return genomes
        genome = [
            self.get_genome_alias(digest=x, fallback=True)
            for x in _make_list_of_str(genome)
        ]
        if strict:
            missing = []
            filtered = []
            for g in genome:
                if g in genomes:
                    filtered.append(g)
                else:
                    missing.append(g)
            if missing:
                _LOGGER.warning(f"Genomes do not include: {', '.join(missing)}")
            return None if not filtered else filtered
        return genomes if not all(x in genomes for x in genome) else genome


def upgrade_config(
    target_version,
    filepath,
    force=False,
    get_json_url=lambda s, i: s + _get_server_endpoints_mapping(s)[i],
    link_fun=lambda s, t: os.symlink(s, t),
):
    """
    Upgrade the config to a selected target version.

    Convert the config file to target_version format, update file structure
    inside genome_folder. Drop genomes for which genome_digest is not available
    on any of the servers and do not have a fasta asset locally.

    :param str target_version: the version updated to
    :param str filepath: path to config file
    :param bool force: whether the upgrade should be confirmed upfront
    :param function(str, str) -> str get_json_url: how to build URL from
            genome server URL base, genome, and asset
    :param callable link_fun: function to use to link files, e.g os.symlink or os.link
    """
    # init rgc obj with provided config
    current_version = yacman.YacAttMap(filepath=filepath)[CFG_VERSION_KEY]

    if current_version == 0.3:
        from .refgenconf_v03 import _RefGenConfV03 as OldRefGenConf

        rgc = OldRefGenConf(filepath=filepath, writable=True)

        if target_version == "0.4":
            from .helpers import alter_file_tree_03_04 as alter_file_tree
            from .helpers import format_config_03_04 as format_config
    else:
        raise NotImplementedError(
            f"Did not upgrade. Upgrade from v{current_version} config is not "
            f"implemented."
        )

    if target_version not in CFG_UPGRADE[str(rgc[CFG_VERSION_KEY])]:
        raise NotImplementedError(
            f"Did not upgrade. Can't upgrade to the requested target "
            f"version ({target_version}). Available target versions for "
            f"{str(rgc[CFG_VERSION_KEY])} are "
            f"{CFG_UPGRADE[str(rgc[CFG_VERSION_KEY])]}"
        )

    # prompt the user
    url = "http://refgenie.databio.org/en/latest/upgrade_config/"
    if not force and not query_yes_no(
        f"Upgrading config to v{target_version}. Current genome identifiers"
        f" will be replaced with sequence-derived digests and contents of "
        f"'{rgc[CFG_FOLDER_KEY]}' will be moved to '{DATA_DIR}' and "
        f"'{ALIAS_DIR}' directories. For more info visit: {url}. Would you "
        f"like to proceed?"
    ):
        _LOGGER.info("Action aborted by the user.")
        return False

    # test server(s) and prompt
    cnt = 0
    incompat_servers = []
    for server in rgc[CFG_SERVERS_KEY]:
        cnt += 1
        try:
            get_json_url(server, API_VERSION + API_ID_ASSETS)
        except (KeyError, ConnectionError, DownloadJsonError):
            incompat_servers.append(server)
    if incompat_servers:
        _LOGGER.info(
            f"The following refgenieserver instances are not "
            f"compatible or do not exist: {incompat_servers}"
        )

    # check digest availability
    missing_digest = []
    for genome, genome_v in rgc[CFG_GENOMES_KEY].items():
        try:
            tag = rgc.get_default_tag(genome, "fasta")
            asset_path = rgc.seek(genome, "fasta", tag, "fasta")
            if not os.path.exists(asset_path):
                raise FileNotFoundError
        except (MissingAssetError, FileNotFoundError):
            cnt = 0
            servers = rgc[CFG_SERVERS_KEY]
            for server in servers:
                cnt += 1
                try:
                    get_json_url(s=server, i=API_ID_ALIAS_DIGEST).format(alias=genome)
                    break
                except (KeyError, ConnectionError, DownloadJsonError) as e:
                    if cnt == len(servers):
                        missing_digest.append(genome)
                    continue

    if (
        not force
        and missing_digest
        and not query_yes_no(
            f"The following genomes will be lost due to the lack of local fasta "
            f"assets and remote genome digests: {', '.join(missing_digest)}. "
            f"Would you like to proceed?"
        )
    ):
        _LOGGER.info("Action aborted by the user.")
        return False

    # reformat config file
    format_config(rgc, get_json_url=get_json_url)
    # alter genome_folder structure
    alter_file_tree(rgc, link_fun=link_fun)
    # change the config_version
    rgc[CFG_VERSION_KEY] = target_version
    # write over the config file
    rgc.write()
    return True


def _download_url_progress(url, output_path, name, params=None):
    """
    Download asset at given URL to given filepath, show progress along the way.

    :param str url: server API endpoint
    :param str output_path: path to file to save download
    :param str name: name to display in front of the progress bar
    :param dict params: query parameters to be added to the request
    """

    class _HookProgress(Progress):
        """
        Internal class to connect progress bar with URL retrieval context manager
        """

        @staticmethod
        def rep_hook(count, blockSize, totalSize):
            """
            Needs to take three arguments in this order
            """
            progress.update(task_id, advance=blockSize)

    def _get_content_len(x):
        """
        Get length of remote content
        """
        f = urlopen(x)
        content_len = f.info().get("Content-length")
        f.close()
        return int(content_len)

    progress = _HookProgress(
        TextColumn("{task.fields[n]}", justify="right"),
        BarColumn(bar_width=None),
        "[magenta]{task.percentage:>3.1f}%",
        "•",
        _DownloadColumn(),
        "•",
        _TransferSpeedColumn(),
        "•",
        _TimeRemainingColumn(),
    )

    url = url if params is None else url + "?{}".format(urlencode(params))
    task_id = progress.add_task("download", n=name, total=_get_content_len(url))
    with progress as p:
        urlretrieve(url, filename=output_path, reporthook=p.rep_hook)


def _genome_asset_path(
    genomes, gname, aname, tname, seek_key, enclosing_dir, no_tag=False
):
    """
    Retrieve the raw path value for a particular asset for a particular genome.

    :param Mapping[str, Mapping[str, Mapping[str, object]]] genomes: nested
        collection of key-value pairs, keyed at top level on genome ID, then by
        asset name, then by asset attribute
    :param str gname: top level key to query -- genome ID, e.g. mm10
    :param str aname: second-level key to query -- asset name, e.g. fasta
    :param str tname: third-level key to query -- tag name, e.g. default
    :param str seek_key: fourth-level key to query -- tag name, e.g. chrom_sizes
    :param bool enclosing_dir: whether a path to the entire enclosing directory should
        be returned, e.g. for a fasta asset that has 3 seek_keys pointing to 3 files
        in an asset dir, that asset dir is returned
    :return str: raw path value for a particular asset for a particular genome
    :raise MissingGenomeError: if the given key-value pair collection does not
        contain as a top-level key the given genome ID
    :raise MissingAssetError: if the given key-value pair collection does
        contain the given genome ID, but that key's mapping doesn't contain
        the given asset name as a key
    :raise GenomeConfigFormatError: if it's discovered during the query that
        the structure of the given genomes mapping suggests that it was
        parsed from an improperly formatted/structured genome config file.
    """
    _assert_gat_exists(genomes, gname, aname, tname)
    asset_tag_data = genomes[gname][CFG_ASSETS_KEY][aname][CFG_ASSET_TAGS_KEY][tname]
    if enclosing_dir:
        if no_tag:
            return asset_tag_data[CFG_ASSET_PATH_KEY]
        return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname)
    if seek_key is None:
        if aname in asset_tag_data[CFG_SEEK_KEYS_KEY]:
            seek_key = aname
        else:
            if no_tag:
                return asset_tag_data[CFG_ASSET_PATH_KEY]
            return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname)
    try:
        seek_key_value = asset_tag_data[CFG_SEEK_KEYS_KEY][seek_key]
    except KeyError:
        raise MissingSeekKeyError(
            f"genome/asset:tag bundle '{gname}/{aname}:{tname}' exists, but "
            f"seek_key '{seek_key}' is missing"
        )
    else:
        if no_tag:
            return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], seek_key_value)
        return os.path.join(asset_tag_data[CFG_ASSET_PATH_KEY], tname, seek_key_value)


def _assert_gat_exists(genomes, gname, aname=None, tname=None, allow_incomplete=False):
    """
    Make sure the genome/asset:tag combination exists in the provided mapping and has
     any seek keys defined. Seek keys are required for the asset completeness.

    :param Mapping[str, Mapping[str, Mapping[str, object]]] genomes: nested
        collection of key-value pairs, keyed at top level on genome ID, then by
        asset name, then by asset attribute
    :param str gname: top level key to query -- genome ID, e.g. mm10
    :param str aname: second-level key to query -- asset name, e.g. fasta
    :param str tname: third-level key to query -- tag name, e.g. default
    :raise MissingGenomeError: if the given key-value pair collection does not
        contain as a top-level key the given genome ID
    :raise MissingAssetError: if the given key-value pair collection does
        contain the given genome ID, but that key's mapping doesn't contain
        the given asset name as a key
    :raise GenomeConfigFormatError: if it's discovered during the query that
        the structure of the given genomes mapping suggests that it was
        parsed from an improperly formatted/structured genome config file.
    """
    _LOGGER.debug(f"checking existence of: {gname}/{aname}:{tname}")
    try:
        genome = genomes[gname]
    except KeyError:
        raise MissingGenomeError(f"Your genomes do not include '{gname}'")
    if aname is not None:
        try:
            asset_data = genome[CFG_ASSETS_KEY][aname]
        except KeyError:
            raise MissingAssetError(
                f"Genome '{gname}' exists, but asset '{aname}' is missing"
            )
        except TypeError:
            _raise_not_mapping(asset_data, "Asset section ")
        if tname is not None:
            try:
                tag_data = asset_data[CFG_ASSET_TAGS_KEY][tname]
            except KeyError:
                raise MissingTagError(
                    f"genome/asset bundle '{gname}/{aname}' exists, "
                    f"but tag '{tname}' is missing"
                )
            except TypeError:
                _raise_not_mapping(asset_data, "Asset section ")
            try:
                tag_data[CFG_SEEK_KEYS_KEY]
            except KeyError:
                if not allow_incomplete:
                    raise MissingSeekKeyError(
                        f"Asset incomplete. No seek keys are defined for "
                        f"'{gname}/{aname}:{tname}'. Build or pull the asset again."
                    )


def _is_large_archive(size, cutoff=10):
    """
    Determines if the file is large based on a string formatted as follows: 15.4GB

    :param str size:  size string
    :return bool: the decision
    """

    def _str2float(x):
        """
        Remove any letters from the file size string and cast the remainder to float
        """
        return float("".join(c for c in x if c in "0123456789."))

    _LOGGER.debug(f"Checking archive size: '{size}'")
    if size.endswith("MB"):
        # convert to gigs
        size = "{0:f}GB".format(_str2float(size) / 1000)
    if size.endswith("KB"):
        # convert to gigs
        size = "{0:f}GB".format(_str2float(size) / 1000 ** 2)
    return size.endswith("TB") or (size.endswith("GB") and _str2float(size) > cutoff)


def _make_genome_assets_line(
    gen,
    assets,
    offset_text="  ",
    genome_assets_delim="/ ",
    asset_sep=", ",
    order=None,
    asset_tag_delim=":",
    rjust=20,
):
    """
    Build a line of text for display of assets by genome

    :param str gen: reference assembly ID, e.g. hg38
    :param Iterable[str] assets: collection of asset names for the given genome
    :param str offset_text: prefix for the line, e.g. a kind of whitespace
    :param str genome_assets_delim: delimiter between a genome ID and text
        showing names of assets for that genome
    :param str asset_sep: delimiter between asset names
    :param function(str) -> object order: how to key asset names for sort
    :return str: text representation of a single assembly's name and assets
    """
    tagged_assets = asset_sep.join(
        sorted(_make_asset_tags_product(assets, asset_tag_delim), key=order)
    )
    return "{}{}{}{}".format(
        gen.rjust(rjust), genome_assets_delim, offset_text, tagged_assets
    )


def _make_asset_tags_product(assets, asset_tag_delim=":", asset_sk_delim="."):
    """
    Make a product of assets and tags available in the provided mapping

    :param Mapping assets: the assets for a selected genome
    :param str asset_tag_delim: how to represent the asset-tag link
    :param str asset_sk_delim: how to represent the asset-seek_key link
    :return list: list representation of tagged assets
    """
    tagged_assets = []
    for aname, asset in assets.items():
        for tname, tag in asset[CFG_ASSET_TAGS_KEY].items():
            sk_assets = []
            seek_keys = get_tag_seek_keys(tag)
            # proceed only if asset is 'complete' -- has seek_keys
            if seek_keys is not None:
                # add seek_keys if exist and different from the asset name,
                # otherwise just the asset name
                sk_assets.extend(
                    [
                        asset_sk_delim.join([aname, sk]) if sk != aname else aname
                        for sk in seek_keys
                    ]
                )
            # add tags to the asset.seek_key list
            tagged_assets.extend(
                [asset_tag_delim.join(i) for i in itertools.product(sk_assets, [tname])]
            )
    return tagged_assets


def _check_insert_data(obj, datatype, name):
    """Checks validity of an object"""
    if obj is None:
        return False
    if not isinstance(obj, datatype):
        raise TypeError(f"{name} must be {datatype.__name__}; got {type(obj).__name__}")
    return True


def _make_list_of_str(arg):
    """
    Convert a str to list of str or ensure a list is a list of str

    :param list[str] | str arg: string or a list of strings to listify
    :return list: list of strings
    :raise TypeError: if a fault argument was provided
    """

    def _raise_faulty_arg():
        raise TypeError(
            f"Provided argument has to be a list[str] or a str, "
            f"got '{arg.__class__.__name__}'"
        )

    if isinstance(arg, str):
        return [arg]
    elif isinstance(arg, list):
        if not all(isinstance(i, str) for i in arg):
            _raise_faulty_arg()
        else:
            return arg
    else:
        _raise_faulty_arg()


def _extend_unique(l1, l2):
    """
    Extend a list with no duplicates

    :param list l1: original list
    :param list l2: list with items to add
    :return list: an extended list
    """
    return l1 + list(set(l2) - set(l1))


def get_asset_tags(asset):
    """
    Return a list of asset tags.

    These need an accession function since under the tag name key there are not only
     tag names, but also the default tag pointer

    :param Mapping asset: a single asset part of the RefGenConf
    :return list: asset tags
    """
    return [t for t in asset[CFG_ASSET_TAGS_KEY]]


def get_tag_seek_keys(tag):
    """
    Return a list of tag seek keys.

    :param Mapping tag: a single tag part of the RefGenConf
    :return list: tag seek keys
    """
    return [s for s in tag[CFG_SEEK_KEYS_KEY]] if CFG_SEEK_KEYS_KEY in tag else None


def construct_request_url(server_url, operation_id, api_prefix=API_VERSION):
    """
    Create a request URL based on a openAPI description

    :param str server_url: server URL
    :param str operation_id: the operationId of the endpoint
    :param str api_prefix: a string to prepend to the operation id
    :return str: a complete URL for the request
    """
    exception_str = f"'{server_url}' is not a compatible refgenieserver instance. "
    try:
        return (
            server_url
            + _get_server_endpoints_mapping(server_url)[api_prefix + operation_id]
        )
    except MissingSchema as e:
        _LOGGER.error(
            exception_str + f"Could not fetch OpenAPI schema: {server_url}/openapi.json"
        )
    except KeyError as e:
        _LOGGER.error(
            exception_str + f"Could not determine API endpoint defined by ID: {e}"
        )


def _get_server_endpoints_mapping(url):
    """
    Establishes the API with the server using operationId field in the openAPI
    JSON description

    :param str url: server URL
    :return dict: endpoints mapped by their operationIds
    """
    json = send_data_request(url + "/openapi.json")
    return map_paths_by_id(
        asciify_json_dict(json) if sys.version_info[0] == 2 else json
    )


def map_paths_by_id(json_dict):
    # check the required input dict characteristics to construct the mapping
    if (
        "openapi" not in json_dict
        or not isinstance(json_dict["openapi"], str)
        or "paths" not in json_dict
        or not isinstance(json_dict["paths"], dict)
    ):
        raise ValueError(
            "The provided mapping is not a valid representation of a "
            "JSON openAPI description"
        )
    return {
        values["get"]["operationId"]: endpoint
        for endpoint, values in json_dict["paths"].items()
    }


def _remove(path):
    """
    remove asset if it is a dir or a file

    :param str path: path to the entity to remove, either a file or a dir
    :return str: removed path
    """
    from shutil import rmtree

    if os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        rmtree(path)
    else:
        raise ValueError(f"path '{path}' is neither a file nor a dir.")
    return path


def _entity_dir_removal_log(directory, entity_class, asset_dict, removed_entities):
    """
    Message and save removed entity data

    :param str directory: removed dir
    :param str entity_class: class of the entity
    :param dict asset_dict: selected genome/asset:tag combination
    :param list removed_entities: list of the removed entities to append to
    """
    subclass = "asset" if entity_class == "genome" else "tag"
    if os.path.basename(directory) == asset_dict[entity_class]:
        _LOGGER.info(
            "Last {sub} for {ec} '{en}' has been removed, "
            "removing {ec} directory".format(
                sub=subclass, ec=entity_class, en=asset_dict[entity_class]
            )
        )
        removed_entities.append(_remove(directory))
    else:
        _LOGGER.debug(
            f"Didn't remove '{directory}' since it does not match the {entity_class} "
            f"name: {asset_dict[entity_class]}"
        )


def _safe_setdef(mapping, attr, val):
    """
    Set default value for a mapping, but catch errors caused by the mapping to
    be updated being an object of incorrect type. Raise an informative error.

    :param Mapping mapping: mapping to update
    :param str attr: attribute to update
    :param val: value to assign as the default
    :raise GenomeConfigFormatError: if mapping is of incorrect class
    :return Mapping: updated mapping
    """
    try:
        mapping.setdefault(attr, val)
    except (TypeError, AttributeError):
        _raise_not_mapping(mapping, f"Cannot update; Section '{attr}' ")
    return mapping


def _raise_not_mapping(mapping, prefix=""):
    raise GenomeConfigFormatError(
        prefix + f"is not a mapping but '{type(mapping).__name__}'. "
        f"This is usually a result of a previous error"
    )


def _populate_refgenie_registry_path(rgc, glob, seek_method_name, remote_class=None):
    """
    Populates refgenie references from refgenie://genome/asset:tag registry paths

    :param dict | str | list glob: String which may contain refgenie registry paths as
        values; or a dict, for which values may contain refgenie registry
        paths. Dict include nested dicts.
    :param bool seek_method_name: a RefGenConf method name to use to seek for the
        strings to replace matched refgenie registry paths, e.g. 'seek' or 'seekr'
    :param str remote_class: remote data provider class. Used only in remote=True
    :return dict | str | list: modified input dict with refgenie paths populated
    """
    p = re.compile("refgenie://([A-Za-z0-9_/\.\:]+)?")
    partial_args = dict()

    # prepare partial function based on operation mode
    if remote_class is not None:
        partial_args.update(dict(remote_class=remote_class))
    _pop = partial(
        rgc.populate if seek_method_name == "seek" else rgc.populater,
        **partial_args,
    )

    if isinstance(glob, str):
        it = re.finditer(p, glob)
        for m in it:
            reg_path = m.group()
            rgpkg = prp(reg_path)
            if not rgpkg:
                _LOGGER.info(
                    f"Can't convert non-conforming refgenie registry path:"
                    f" {reg_path}"
                )
                return glob
            args = dict(
                genome_name=rgpkg["namespace"],
                asset_name=rgpkg["item"],
                tag_name=rgpkg["tag"],
                seek_key=rgpkg["subitem"],
            )
            if remote_class is not None:
                args.update(dict(remote_class=remote_class))
            rgpath = getattr(rgc, seek_method_name)(**args)
            if rgpath is None:
                _LOGGER.warning(f"'{reg_path}' refgenie registry path not populated.")
                continue
            glob = re.sub(reg_path, rgpath, glob)
        return glob
    elif isinstance(glob, dict):
        for k, v in glob.items():
            if k.startswith("_"):
                continue
            if k.startswith("sources"):
                continue  # derived attribute sources
            glob[k] = _pop(v)
        return glob
    elif isinstance(glob, list):
        return [_pop(v) for v in glob]
    elif isinstance(glob, AttMap):
        return AttMap(_pop(glob.to_dict()))
    else:
        otype = type(glob)
        _LOGGER.debug(
            f"Refgenie can only populate str, list, or dict objects. Got {otype}"
        )
        _LOGGER.debug(glob)
        return glob
