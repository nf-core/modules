import json
import os
import sys
from collections import OrderedDict
from functools import partial

import logmuse
from refgenconf import (
    DownloadJsonError,
    MissingAssetError,
    MissingGenomeError,
    RefGenConf,
)
from refgenconf import __version__ as rgc_version
from refgenconf import select_genome_config, upgrade_config
from refgenconf.helpers import block_iter_repr
from requests.exceptions import MissingSchema
from rich.console import Console
from ubiquerg import query_yes_no

from ._version import __version__
from .argparser import build_argparser
from .asset_build_packages import *
from .const import *
from .exceptions import *
from .helpers import _raise_missing_recipe_error, _single_folder_writeable
from .refgenie import (
    _skip_lock,
    parse_registry_path,
    refgenie_build,
    refgenie_build_reduce,
)


def main():
    """Primary workflow"""
    parser = logmuse.add_logging_options(build_argparser())
    args, _ = parser.parse_known_args()
    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args, make_root=True)
    _LOGGER.debug(f"versions: refgenie {__version__} | refgenconf {rgc_version}")
    _LOGGER.debug(f"Args: {args}")

    if not args.command:
        parser.print_help()
        _LOGGER.error("No command given")
        sys.exit(1)

    if args.command == ALIAS_CMD and not args.subcommand:
        parser.print_help()
        _LOGGER.error("No alias subcommand command given")
        sys.exit(1)

    if (
        args.command == BUILD_CMD
        and args.asset_registry_paths is None
        and not args.reduce
    ):
        parser.error("You must provide an asset-registry-path")
        sys.exit(1)

    gencfg = select_genome_config(
        filename=args.genome_config,
        check_exist=not args.command == INIT_CMD,
        on_missing=lambda fp: fp,
        strict_env=True,
    )
    if gencfg is None and args.command not in [
        GET_REMOTE_ASSET_CMD,
        LIST_REMOTE_CMD,
        POPULATE_REMOTE_CMD,
    ]:
        raise MissingGenomeConfigError(args.genome_config)
    _LOGGER.debug("Determined genome config: {}".format(gencfg))

    skip_read_lock = True if gencfg is None else _skip_lock(args.skip_read_lock, gencfg)
    # From user input we want to construct a list of asset dicts, where each
    # asset has a genome name, asset name, and tag
    if "asset_registry_paths" in args and args.asset_registry_paths:
        _LOGGER.debug("Found registry_path: {}".format(args.asset_registry_paths))
        asset_list = [parse_registry_path(x) for x in args.asset_registry_paths]
        # [{"protocol": 'pname', "genome": 'gname', "asset", 'aname', "seek_key", 'sname', "tag": 'tname'}, ...]
        for a in asset_list:
            # every asset must have a genome, either provided via registry path
            # or the args.genome arg.
            if not a["genome"]:
                if args.genome:
                    a["genome"] = args.genome
                else:
                    _LOGGER.error(
                        "Provided asset registry path ({}/{}:{}) is invalid. See help for usage reference.".format(
                            a["genome"], a["asset"], a["tag"]
                        )
                    )
                    sys.exit(1)
            else:
                if args.genome and args.genome != a["genome"]:
                    _LOGGER.warn(
                        "Two different genomes specified for asset '{}'.".format(
                            a["asset"]
                        )
                    )

    if args.command == INIT_CMD:
        _LOGGER.debug("Initializing refgenie genome configuration")
        entries = OrderedDict(
            {
                CFG_VERSION_KEY: REQ_CFG_VERSION,
                CFG_FOLDER_KEY: os.path.dirname(os.path.abspath(gencfg)),
                CFG_SERVERS_KEY: args.genome_server or [DEFAULT_SERVER],
                CFG_GENOMES_KEY: None,
            }
        )
        if args.settings_json:
            if os.path.isfile(args.settings_json):
                with open(args.settings_json, "r") as json_file:
                    data = json.load(json_file)
                entries.update(data)
            else:
                raise FileNotFoundError(
                    "JSON file with config init settings does not exist: {}".format(
                        args.settings_json
                    )
                )
        if args.genome_folder:
            entries.update({CFG_FOLDER_KEY: args.genome_folder})
        if args.genome_archive_folder:
            entries.update({CFG_ARCHIVE_KEY: args.genome_archive_folder})
        if args.genome_archive_config:
            entries.update({CFG_ARCHIVE_CONFIG_KEY: args.genome_archive_config})
        _LOGGER.debug("initializing with entries: {}".format(entries))
        rgc = RefGenConf(entries=entries, skip_read_lock=skip_read_lock)
        rgc.initialize_config_file(os.path.abspath(gencfg))

    elif args.command == BUILD_CMD:
        if args.reduce:
            refgenie_build_reduce(
                gencfg=gencfg,
                preserve_map_configs=args.preserve_map_configs,
            )
            sys.exit(0)
        if not all([x["genome"] == asset_list[0]["genome"] for x in asset_list]):
            _LOGGER.error("Build can only build assets for one genome")
            sys.exit(1)
        recipe_name = None
        if args.recipe:
            if len(asset_list) > 1:
                _LOGGER.error("Recipes cannot be specified for multi-asset builds")
                sys.exit(1)
            recipe_name = args.recipe
        if args.requirements:
            for a in asset_list:
                recipe = recipe_name or a["asset"]
                if recipe not in asset_build_packages.keys():
                    _raise_missing_recipe_error(recipe)
                _LOGGER.info("'{}' recipe requirements: ".format(recipe))
                _make_asset_build_reqs(recipe)
            sys.exit(0)

        ret = refgenie_build(
            gencfg, asset_list[0]["genome"], asset_list, recipe_name, args
        )
        if not ret:
            sys.exit(1)
        else:
            sys.exit(0)

    elif args.command == GET_ASSET_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        check = args.check_exists if args.check_exists else None
        for a in asset_list:
            _LOGGER.debug(
                "getting asset: '{}/{}.{}:{}'".format(
                    a["genome"], a["asset"], a["seek_key"], a["tag"]
                )
            )
            print(
                rgc.seek(
                    a["genome"],
                    a["asset"],
                    a["tag"],
                    a["seek_key"],
                    strict_exists=check,
                )
            )
        return

    elif args.command == GET_REMOTE_ASSET_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        if args.genome_server is not None:
            rgc.subscribe(
                urls=args.genome_server, reset=not args.append_server, no_write=True
            )
        for a in asset_list:
            _LOGGER.debug(
                "getting remote asset path: '{}/{}.{}:{}'".format(
                    a["genome"], a["asset"], a["seek_key"], a["tag"]
                )
            )
            print(
                rgc.seekr(
                    a["genome"],
                    a["asset"],
                    a["tag"],
                    a["seek_key"],
                    args.remote_class,
                )
            )
        return

    elif args.command == INSERT_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        if len(asset_list) > 1:
            raise NotImplementedError("Can only add 1 asset at a time")
        else:
            sk = args.seek_keys
            if sk:
                sk = json.loads(args.seek_keys)
            rgc.add(
                path=args.path,
                genome=asset_list[0]["genome"],
                asset=asset_list[0]["asset"],
                tag=asset_list[0]["tag"],
                seek_keys=sk,
                force=args.force,
            )

    elif args.command == PULL_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)

        # existing assets overwriting
        if args.no_overwrite:
            force = False
        elif args.force_overwrite:
            force = True
        else:
            force = None
        # large archive pulling
        if args.no_large:
            force_large = False
        elif args.pull_large:
            force_large = True
        else:
            force_large = None
        # batch mode takes precedence over other choices
        if args.batch:
            force_large = True
            force = False

        outdir = rgc.data_dir
        if not os.path.exists(outdir):
            raise MissingFolderError(outdir)
        if not perm_check_x(outdir):
            return
        if not _single_folder_writeable(outdir):
            _LOGGER.error("Insufficient permissions to write to: {}".format(outdir))
            return

        for a in asset_list:
            rgc.pull(
                a["genome"],
                a["asset"],
                a["tag"],
                force=force,
                force_large=force_large,
                size_cutoff=args.size_cutoff,
            )

    elif args.command in [LIST_LOCAL_CMD, LIST_REMOTE_CMD]:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        console = Console()
        if args.command == LIST_REMOTE_CMD:
            if args.genome_server is not None:
                rgc.subscribe(
                    urls=args.genome_server, reset=not args.append_server, no_write=True
                )
            num_servers = 0
            bad_servers = []
            for server_url in rgc[CFG_SERVERS_KEY]:
                num_servers += 1
                try:
                    table = rgc.get_asset_table(
                        genomes=args.genome, server_url=server_url
                    )
                except (DownloadJsonError, ConnectionError, MissingSchema):
                    bad_servers.append(server_url)
                    continue
                else:
                    console.print(table)
            if num_servers >= len(rgc[CFG_SERVERS_KEY]) and bad_servers:
                _LOGGER.error(
                    "Could not list assets from the following servers: {}".format(
                        bad_servers
                    )
                )
        else:
            if args.recipes:
                print(", ".join(sorted(list(asset_build_packages.keys()))))
            else:
                console.print(rgc.get_asset_table(genomes=args.genome))

    elif args.command == GETSEQ_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        print(rgc.getseq(args.genome, args.locus))

    elif args.command == REMOVE_CMD:
        force = args.force
        rgc = RefGenConf(filepath=gencfg, skip_read_lock=skip_read_lock)
        for a in asset_list:
            a["tag"] = a["tag"] or rgc.get_default_tag(
                a["genome"], a["asset"], use_existing=False
            )
            _LOGGER.debug("Determined tag for removal: {}".format(a["tag"]))
            if a["seek_key"] is not None:
                raise NotImplementedError("You can't remove a specific seek_key.")
            gat = {"genome": a["genome"], "asset": a["asset"], "tag": a["tag"]}
            try:
                if not rgc.is_asset_complete(**gat):
                    with rgc as r:
                        r.cfg_remove_assets(**gat)
                    _LOGGER.info(
                        "Removed an incomplete asset "
                        "'{genome}/{asset}:{tag}'".format(*gat)
                    )
                    return
            except (KeyError, MissingAssetError, MissingGenomeError):
                _LOGGER.info(
                    "Asset '{genome}/{asset}:{tag}' does not exist".format(**gat)
                )
                return
        if len(asset_list) > 1:
            if not query_yes_no(
                "Are you sure you want to remove {} assets?".format(len(asset_list))
            ):
                _LOGGER.info("Action aborted by the user")
                return
            force = True
        for a in asset_list:
            rgc.remove(genome=a["genome"], asset=a["asset"], tag=a["tag"], force=force)

    elif args.command == TAG_CMD:
        rgc = RefGenConf(filepath=gencfg, skip_read_lock=skip_read_lock)
        if len(asset_list) > 1:
            raise NotImplementedError("Can only tag 1 asset at a time")
        if args.default:
            # set the default tag and exit
            with rgc as r:
                r.set_default_pointer(a["genome"], a["asset"], a["tag"], True)
            sys.exit(0)
        rgc.tag(a["genome"], a["asset"], a["tag"], args.tag, force=args.force)

    elif args.command == ID_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        if len(asset_list) == 1:
            g, a = asset_list[0]["genome"], asset_list[0]["asset"]
            t = asset_list[0]["tag"] or rgc.get_default_tag(g, a)
            print(rgc.id(g, a, t))
            return
        for asset in asset_list:
            g, a = asset["genome"], asset["asset"]
            t = asset["tag"] or rgc.get_default_tag(g, a)
            print("{}/{}:{},".format(g, a, t) + rgc.id(g, a, t))
        return
    elif args.command == SUBSCRIBE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        rgc.subscribe(urls=args.genome_server, reset=args.reset)
        return
    elif args.command == UNSUBSCRIBE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        rgc.unsubscribe(urls=args.genome_server)
        return
    elif args.command == ALIAS_CMD:
        rgc = RefGenConf(filepath=gencfg, skip_read_lock=skip_read_lock)
        if args.subcommand == ALIAS_GET_CMD:
            if args.aliases is not None:
                for a in args.aliases:
                    print(rgc.get_genome_alias_digest(alias=a))
                return
            console = Console()
            console.print(rgc.genome_aliases_table)

        if args.subcommand == ALIAS_SET_CMD:
            rgc.set_genome_alias(
                digest=args.digest,
                genome=args.aliases,
                reset_digest=args.reset,
                create_genome=args.force,
            )
            return
        elif args.subcommand == ALIAS_REMOVE_CMD:
            rgc.remove_genome_aliases(digest=args.digest, aliases=args.aliases)
            return

    elif args.command == COMPARE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        res = rgc.compare(
            args.genome1[0], args.genome2[0], explain=not args.no_explanation
        )
        if args.no_explanation:
            print(res)
        if args.flag_meanings:
            from refgenconf.seqcol import FLAGS
            from rich.table import Table

            _LOGGER.info("\n")
            codes = sorted(FLAGS.keys())
            table = Table(title="Compatibility flags")
            table.add_column("Code")
            table.add_column("Indication")
            for code in codes:
                table.add_row(str(code), FLAGS[code])
            console = Console()
            console.print(table)

    elif args.command == UPGRADE_CMD:
        upgrade_config(
            target_version=args.target_version, filepath=gencfg, force=args.force
        )

    elif args.command == POPULATE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        process_populate(pop_fun=rgc.populate, file_path=args.file)

    elif args.command == POPULATE_REMOTE_CMD:
        rgc = RefGenConf(filepath=gencfg, writable=False, skip_read_lock=skip_read_lock)
        if args.genome_server is not None:
            rgc.subscribe(
                urls=args.genome_server, reset=not args.append_server, no_write=True
            )
        pop_fun = partial(rgc.populater, remote_class=args.remote_class)
        process_populate(pop_fun=pop_fun, file_path=args.file)


def process_populate(pop_fun, file_path=None):
    """
    Process a populate request (file or stdin) with a custom populator function

    :param callable(dict | str | list) -> dict | str | list pop_fun: a function
        that populates refgenie registry paths in objects
    :param str file_path: path to the file to populate refgenie registry paths in,
        skip for stdin processing
    """
    if file_path is not None:
        _LOGGER.debug(f"Populating file: {file_path}")
        with open(file_path) as fp:
            for line in fp:
                sys.stdout.write(pop_fun(glob=line))
    else:
        for line in sys.stdin:
            if line.rstrip() in ["q", "quit", "exit"]:
                break
            sys.stdout.write(pop_fun(glob=line))


def perm_check_x(file_to_check, message_tag="genome directory"):
    """
    Check X_OK permission on a path, providing according messaging and bool val.

    :param str file_to_check: path to query for permission
    :param str message_tag: context for error message if check fails
    :return bool: os.access(path, X_OK) for the given path
    :raise ValueError: if there's no filepath to check for permission
    """
    if not file_to_check:
        msg = "You must provide a path to {}".format(message_tag)
        _LOGGER.error(msg)
        raise ValueError(msg)
    if not os.access(file_to_check, os.X_OK):
        _LOGGER.error("Insufficient permissions to write to {}: ".format(file_to_check))
        return False
    return True


def _make_asset_build_reqs(asset):
    """
    Prepare requirements and inputs lists and display it

    :params str asset: name of the asset
    """

    def _format_reqs(req_list):
        """

        :param list[dict] req_list:
        :return list[str]:
        """
        templ = "\t{} ({})"
        return [
            templ.format(req[KEY], req[DESC])
            if DEFAULT not in req
            else (templ + "; default: {}").format(req[KEY], req[DESC], req[DEFAULT])
            for req in req_list
        ]

    reqs_list = []
    if asset_build_packages[asset][REQ_FILES]:
        reqs_list.append(
            "- files:\n{}".format(
                "\n".join(_format_reqs(asset_build_packages[asset][REQ_FILES]))
            )
        )
    if asset_build_packages[asset][REQ_ASSETS]:
        reqs_list.append(
            "- assets:\n{}".format(
                "\n".join(_format_reqs(asset_build_packages[asset][REQ_ASSETS]))
            )
        )
    if asset_build_packages[asset][REQ_PARAMS]:
        reqs_list.append(
            "- params:\n{}".format(
                "\n".join(_format_reqs(asset_build_packages[asset][REQ_PARAMS]))
            )
        )
    _LOGGER.info("\n".join(reqs_list))
