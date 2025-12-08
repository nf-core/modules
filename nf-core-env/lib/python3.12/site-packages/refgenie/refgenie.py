import csv
import json
import os
import re
import signal
import sys
from glob import glob
from logging import getLogger

import pypiper
from refgenconf import RefGenConf, get_dir_digest
from refgenconf.exceptions import (
    MissingAssetError,
    MissingGenomeError,
    MissingSeekKeyError,
    MissingTagError,
)
from refgenconf.helpers import block_iter_repr
from rich.progress import track
from ubiquerg import parse_registry_path as prp
from ubiquerg.system import is_writable
from yacman import UndefinedAliasError

from .asset_build_packages import *
from .const import *
from .helpers import (
    _parse_user_build_input,
    _raise_missing_recipe_error,
    _skip_lock,
    _writeable,
    make_sure_path_exists,
)

_LOGGER = getLogger(PKG_NAME)


def parse_registry_path(path):
    return prp(
        path,
        defaults=[
            ("protocol", None),
            ("genome", None),
            ("asset", None),
            ("seek_key", None),
            ("tag", None),
        ],
    )


def get_asset_vars(
    genome,
    asset_key,
    tag,
    outfolder,
    specific_args=None,
    specific_params=None,
    **kwargs,
):
    """
    Gives a dict with variables used to populate an asset path.
    """
    asset_outfolder = os.path.join(outfolder, asset_key, tag)
    asset_vars = {
        "genome": genome,
        "asset": asset_key,
        "tag": tag,
        "asset_outfolder": asset_outfolder,
    }
    if specific_args:
        asset_vars.update(specific_args)
    if specific_params:
        asset_vars.update(specific_params)
    asset_vars.update(**kwargs)
    return asset_vars


def refgenie_initg(rgc, genome, content_checksums):
    """
    Initializing a genome means adding `collection_checksum` attributes in the
    genome config file. This should perhaps be a function in refgenconf, but not
    a CLI-hook. Also adds `content_checksums` tsv file (should be a recipe cmd?).

    This function updates the provided RefGenConf object with the
    genome(collection)-level checksum and saves the individual checksums to a
    TSV file in the fasta asset directory.

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param str genome: name of the genome
    :param dict content_checksums: checksums of individual content_checksums, e.g. chromosomes
    """
    genome_dir = os.path.join(rgc.data_dir, genome)
    if is_writable(genome_dir):
        output_file = os.path.join(genome_dir, "{}_sequence_digests.tsv".format(genome))
        with open(output_file, "w") as contents_file:
            wr = csv.writer(contents_file, delimiter="\t")
            for key, val in content_checksums.items():
                wr.writerow([key, val])
        _LOGGER.debug("sequence digests saved to: {}".format(output_file))
    else:
        _LOGGER.warning(
            "Could not save the genome sequence digests. '{}' is not writable".format(
                genome_dir
            )
        )


def refgenie_build_reduce(gencfg, preserve_map_configs=False):
    """
    Asset building process may be split into two tasks: building assets (_Map_ procedure)
    and gathering asset metadata (_Reduce_ procedure).

    This function performs the _Reduce_ procedure:
    finds the genome configuration files produced in the _Map_ step,
    updates the main genome configuration file with their contents and removes them.

    :param str gencfg: an absolute path to the genome configuration file
    :param bool preserve_map_configs: a boolean indicating whether the map configs should be preserved,
        by default they are removed once the contents are integrated into the master genome config.
    :return bool: a boolean indicating whether the master config has been successfully updated
        or None in case there were no map configs found.
    """

    def _map_cfg_match_pattern(data_dir, match_all_str):
        """
        Create a path to the map genome config witb a provided 'match all' character,
        which needs to be different depending on the matchig scenario.

        :param str data_dir: an absolute path to the data directory
        :param str match_all_str: match all character to use
        """
        return os.path.join(
            data_dir,
            *([match_all_str] * 3),
            BUILD_STATS_DIR,
            BUILD_MAP_CFG,
        )

    _LOGGER.info("Running the reduce procedure. No assets will be built.")
    rgc_master = RefGenConf(filepath=gencfg, writable=True)
    regex_pattern = _map_cfg_match_pattern(rgc_master.data_dir, "(\S+)")
    glob_pattern = _map_cfg_match_pattern(rgc_master.data_dir, "*")
    rgc_map_filepaths = glob(glob_pattern, recursive=True)
    if len(rgc_map_filepaths) == 0:
        _LOGGER.info(f"No map configs to reduce")
        return None
    _LOGGER.debug(f"Map configs to reduce: {block_iter_repr(rgc_map_filepaths)}")
    matched_gats = []
    for rgc_map_filepath in track(
        rgc_map_filepaths,
        description=f"Reducing {len(rgc_map_filepaths)} configs",
    ):
        matched_genome, matched_asset, matched_tag = re.match(
            pattern=regex_pattern, string=rgc_map_filepath
        ).groups()
        matched_gat = f"{matched_genome}/{matched_asset}:{matched_tag}"
        map_rgc = RefGenConf(filepath=rgc_map_filepath, writable=False)
        if CFG_GENOMES_KEY not in map_rgc:
            _LOGGER.warning(
                f"'{rgc_map_filepath}' is missing '{CFG_GENOMES_KEY}' key, skipping"
            )
            continue
        # this should be a one element list
        genome_digests = map_rgc[CFG_GENOMES_KEY].keys()
        if len(genome_digests) > 1:
            _LOGGER.warning(
                f"There are {len(genome_digests)} genomes in the map build config while 1 expected, skipping"
            )
            continue
        genome_digest = genome_digests[0]
        alias = map_rgc.get_genome_alias(digest=genome_digest)
        if genome_digest != matched_genome:
            raise Exception(
                f"Genome directory name does not match genome in the map config: {matched_genome} != {genome_digest}"
            )
        asset_data = tag_data = map_rgc[CFG_GENOMES_KEY][matched_genome][
            CFG_ASSETS_KEY
        ][matched_asset]
        tag_data = asset_data[CFG_ASSET_TAGS_KEY][matched_tag]
        default_tag_in_map = asset_data[CFG_ASSET_DEFAULT_TAG_KEY]
        try:
            alias_master = rgc_master.get_genome_alias(digest=genome_digest)
            assert alias == alias_master
        except (UndefinedAliasError, AssertionError):
            # no need to put this in context manager
            # it is already used in the method
            rgc_master.set_genome_alias(
                genome=alias, digest=genome_digest, create_genome=True
            )
        with rgc_master as r:
            if CFG_ASSET_PARENTS_KEY in tag_data:
                for parent in tag_data[CFG_ASSET_PARENTS_KEY]:
                    parsed_parent = parse_registry_path(parent)
                    r.update_relatives_assets(
                        genome=parsed_parent["genome"],
                        asset=parsed_parent["asset"],
                        tag=parsed_parent["tag"],
                        data=[matched_gat],
                        children=True,
                    )

            if CFG_ASSET_CHILDREN_KEY in tag_data:
                for child in tag_data[CFG_ASSET_CHILDREN_KEY]:
                    parsed_child = parse_registry_path(child)
                    r.update_relatives_assets(
                        genome=parsed_child["genome"],
                        asset=parsed_child["asset"],
                        tag=parsed_child["tag"],
                        data=[matched_gat],
                        children=False,
                    )
            r.update_tags(
                genome=matched_genome,
                asset=matched_asset,
                tag=matched_tag,
                data=tag_data,
                force_digest=genome_digest,
            )
            # set a default tag in the master config to the one built in map mode,
            # this will not overwrite an existing tag though
            r.set_default_pointer(
                genome=matched_genome,
                asset=matched_asset,
                tag=default_tag_in_map,
            )
        matched_gats.append(matched_gat)
        if not preserve_map_configs:
            os.remove(rgc_map_filepath)
    _LOGGER.info(f"Added entries for: {block_iter_repr(matched_gats)}")
    return True


def refgenie_build(gencfg, genome, asset_list, recipe_name, args):
    """
    Runs the refgenie build recipe.

    :param str gencfg: path to the genome configuration file
    :param str genome:
    :param list asset_list:
    :param str recipe_name:
    :param argparse.Namespace args: parsed command-line options/arguments
    """

    rgc = RefGenConf(
        filepath=gencfg,
        writable=False,
        skip_read_lock=_skip_lock(args.skip_read_lock, gencfg),
    )
    specified_args = _parse_user_build_input(args.files)
    specified_params = _parse_user_build_input(args.params)

    def _read_json_file(filepath):
        """
        Read a JSON file

        :param str filepath: path to the file to read
        :return dict: read data
        """
        with open(filepath, "r") as f:
            data = json.load(f)
        return data

    if recipe_name and os.path.isfile(recipe_name) and recipe_name.endswith(".json"):
        recipe_name = _read_json_file(filepath=recipe_name)

    def _build_asset(
        genome,
        asset_key,
        tag,
        build_pkg,
        genome_outfolder,
        specific_args,
        specific_params,
        alias,
        **kwargs,
    ):
        """
        Builds assets with pypiper and updates a genome config file.

        This function actually runs the build commands in a given build package,
        and then update the refgenie config file.

        :param str genome: The assembly key; e.g. 'mm10'.
        :param str asset_key: The unique asset identifier; e.g. 'bowtie2_index'
        :param dict build_pkg: A dict (see examples) specifying lists
            of required input_assets, commands to run, and outputs to register as
            assets.
        """
        if args.map:
            # Performing a build map step.
            # The reduce step will need to be performed to get the built
            # asset metadata to the master config file
            genome_alias = rgc.get_genome_alias(digest=genome)
            # create an empty config file in the genome directory
            _LOGGER.info(f"Using new map genome config: {locked_map_gencfg}")
            make_sure_path_exists(os.path.dirname(locked_map_gencfg))
            open(locked_map_gencfg, "a").close()
            # initialize a new RefGenConf.
            # Use the master location for data storage,
            # but change path to the in asset dir location
            rgc_map = RefGenConf(
                entries={"genome_folder": rgc.genome_folder},
                filepath=locked_map_gencfg,
            )
            # set the alias first (if available), based on the master file

            rgc_map.set_genome_alias(
                digest=genome,
                genome=genome_alias,
                create_genome=True,
            )

            # copy the genome of interest section to the new RefGenConf,
            # so that possible dependancies can be satisfied
            rgc_map.update_genomes(
                genome=genome_alias,
                data=rgc[CFG_GENOMES_KEY][genome],
            )

        else:
            rgc_map = rgc

        _LOGGER.info(
            f"Saving outputs to:{block_iter_repr(['content: ' + genome_outfolder, 'logs: ' + build_stats_dir])}"
        )
        if args.docker:
            # Set up some docker stuff
            if args.volumes:
                # TODO: is volumes list defined here?
                volumes = volumes.append(genome_outfolder)
            else:
                volumes = genome_outfolder

        if not _writeable(genome_outfolder):
            _LOGGER.error(
                f"Insufficient permissions to write to output folder: {genome_outfolder}"
            )
            return False, rgc_map

        pm = pypiper.PipelineManager(
            name=PKG_NAME, outfolder=build_stats_dir, args=args
        )
        tk = pypiper.NGSTk(pm=pm)
        if args.docker:
            pm.get_container(build_pkg[CONT], volumes)
        _LOGGER.debug("Asset build package: " + str(build_pkg))
        # create a bundle list to simplify calls below
        gat = [genome, asset_key, tag]
        # collect variables required to populate the command templates
        asset_vars = get_asset_vars(
            genome,
            asset_key,
            tag,
            genome_outfolder,
            specific_args,
            specific_params,
            **kwargs,
        )
        # populate command templates
        # prior to populating, remove any seek_key parts from the keys, since these are not supported by format method
        command_list_populated = [
            x.format(**{k.split(".")[0]: v for k, v in asset_vars.items()})
            for x in build_pkg[CMD_LST]
        ]
        # create output directory
        tk.make_dir(asset_vars["asset_outfolder"])

        target = os.path.join(
            build_stats_dir, TEMPLATE_TARGET.format(genome, asset_key, tag)
        )
        # add target command
        command_list_populated.append("touch {target}".format(target=target))
        _LOGGER.debug(
            "Command populated: '{}'".format(" ".join(command_list_populated))
        )
        try:
            # run build command
            signal.signal(signal.SIGINT, _handle_sigint(gat))
            pm.run(command_list_populated, target, container=pm.container)
        except pypiper.exceptions.SubprocessError:
            _LOGGER.error("asset '{}' build failed".format(asset_key))
            return False, rgc_map
        else:
            # save build recipe to the JSON-formatted file
            recipe_file_name = TEMPLATE_RECIPE_JSON.format(asset_key, tag)
            with open(os.path.join(build_stats_dir, recipe_file_name), "w") as outfile:
                json.dump(build_pkg, outfile)
            # since the assets are always built to a standard dir structure, we
            # can just stitch a path together for asset digest calculation
            asset_dir = os.path.join(rgc_map.data_dir, *gat)
            if not os.path.exists(asset_dir):
                raise OSError(
                    "Could not compute asset digest. Path does not "
                    "exist: {}".format(asset_dir)
                )
            digest = get_dir_digest(asset_dir)
            _LOGGER.info(f"Asset digest: {digest}")
            # add a 'dir' seek_key that points to the asset directory
            build_pkg[ASSETS].update({"dir": "."})
            # add updates to config file
            with rgc_map as r:
                if asset_key == "fasta":
                    r.update_genomes(
                        genome, data={CFG_ALIASES_KEY: [alias]}, force_digest=genome
                    )
                r.update_assets(
                    *gat[0:2],
                    data={CFG_ASSET_DESC_KEY: build_pkg[DESC]},
                    force_digest=genome,
                )
                r.update_tags(
                    *gat,
                    force_digest=genome,
                    data={
                        CFG_ASSET_PATH_KEY: asset_key,
                        CFG_ASSET_CHECKSUM_KEY: digest,
                    },
                )
                r.update_seek_keys(
                    *gat,
                    force_digest=genome,
                    keys={
                        k: v.format(**asset_vars) for k, v in build_pkg[ASSETS].items()
                    },
                )
                r.set_default_pointer(*gat, force_digest=genome)
        pm.stop_pipeline()
        return True, rgc_map

    for a in asset_list:
        asset_key = a["asset"]
        asset_tag = a["tag"] or rgc.get_default_tag(
            genome, a["asset"], use_existing=False
        )
        recipe_name = recipe_name or asset_key

        if isinstance(recipe_name, dict) or (
            isinstance(recipe_name, str) and recipe_name in asset_build_packages.keys()
        ):
            if isinstance(recipe_name, dict):
                _LOGGER.info("Using custom recipe: \n{}".format(recipe_name))
                asset_build_package = _check_recipe(recipe_name)
                recipe_name = asset_build_package["name"]
            else:
                asset_build_package = _check_recipe(asset_build_packages[recipe_name])
            # handle user-requested parents for the required assets
            input_assets = {}
            parent_assets = []
            specified_asset_keys, specified_assets = None, None
            if args.assets is not None:
                parsed_parents_input = _parse_user_build_input(args.assets)
                specified_asset_keys = list(parsed_parents_input.keys())
                specified_assets = list(parsed_parents_input.values())
                _LOGGER.debug(f"Custom assets requested: {args.assets}")
            if not specified_asset_keys and isinstance(args.assets, list):
                _LOGGER.warning(
                    "Specified parent assets format is invalid. Using defaults."
                )
            for req_asset in asset_build_package[REQ_ASSETS]:
                req_asset_data = parse_registry_path(req_asset[KEY])
                # for each req asset see if non-default parents were requested
                if (
                    specified_asset_keys is not None
                    and req_asset_data["asset"] in specified_asset_keys
                ):
                    parent_data = parse_registry_path(
                        specified_assets[
                            specified_asset_keys.index(req_asset_data["asset"])
                        ]
                    )
                    g, a, t, s = (
                        parent_data["genome"],
                        parent_data["asset"],
                        parent_data["tag"]
                        or rgc.get_default_tag(genome, parent_data["asset"]),
                        parent_data["seek_key"],
                    )
                else:  # if no custom parents requested for the req asset, use default one
                    default = parse_registry_path(req_asset[DEFAULT])
                    g, a, t, s = (
                        genome,
                        default["asset"],
                        rgc.get_default_tag(genome, default["asset"]),
                        req_asset_data["seek_key"],
                    )
                try:
                    parent_assets.append(
                        "{}/{}:{}".format(
                            rgc.get_genome_alias_digest(g, fallback=True), a, t
                        )
                    )
                except UndefinedAliasError as e:
                    _LOGGER.warning(f"'{g}' namespace has not been initialized yet")
                    if args.pull_parents:
                        _LOGGER.info(f"Pulling missing parent: {g}/{a}:{t}")
                        ret = rgc.pull(genome=g, asset=a, tag=t)
                        if ret is None or not all(ret):
                            _LOGGER.info(
                                f"Missing parent asset pull requested, but failed: {g}/{a}:{t}. "
                                f"Reason: {str(e)}"
                            )
                            return False
                    else:
                        raise
                try:
                    input_assets[req_asset[KEY]] = _seek(rgc, g, a, t, s)
                except (
                    MissingAssetError,
                    MissingGenomeError,
                    MissingTagError,
                    MissingSeekKeyError,
                ) as e:
                    if args.pull_parents:
                        _LOGGER.info(f"Pulling missing parent: {g}/{a}:{t}")
                        ret = rgc.pull(genome=g, asset=a, tag=t)
                        if ret is None or not all(ret):
                            _LOGGER.info(
                                f"Missing parent asset pull requested, but failed: {g}/{a}:{t}. "
                                f"Reason: {str(e)}"
                            )
                            return False
                    else:
                        raise

            _LOGGER.debug("Using parents: {}".format(", ".join(parent_assets)))
            _LOGGER.debug("Provided files: {}".format(specified_args))
            _LOGGER.debug("Provided parameters: {}".format(specified_params))
            for required_file in asset_build_package[REQ_FILES]:
                if (
                    specified_args is None
                    or required_file[KEY] not in specified_args.keys()
                ):
                    raise ValueError(
                        "Path to the '{x}' input ({desc}) is required, but not provided. "
                        "Specify it with: --files {x}=/path/to/{x}_file".format(
                            x=required_file[KEY], desc=required_file[DESC]
                        )
                    )
            for required_param in asset_build_package[REQ_PARAMS]:
                if specified_params is None:
                    specified_params = {}
                if required_param[KEY] not in specified_params.keys():
                    if required_param[DEFAULT] is None:
                        raise ValueError(
                            "Value for the parameter '{x}' ({desc}) is required, but not provided. "
                            "Specify it with: --params {x}=value".format(
                                x=required_param[KEY], desc=required_param[DESC]
                            )
                        )
                    else:
                        specified_params.update(
                            {required_param[KEY]: required_param[DEFAULT]}
                        )
            _LOGGER.info(
                "Building '{}/{}:{}' using '{}' recipe".format(
                    genome, asset_key, asset_tag, recipe_name
                )
            )
            ori_genome = genome
            if recipe_name == "fasta":
                if (
                    genome in rgc.genomes_list()
                    and "fasta" in rgc.list_assets_by_genome(genome)
                ):
                    pretag = rgc.get_default_tag(genome, "fasta")
                    _LOGGER.warning(
                        "'{g}' genome is already initialized with other fasta asset ({g}/{a}:{t})".format(
                            g=genome, a=asset_key, t=pretag
                        )
                    )
                    genome = rgc.get_genome_alias_digest(alias=genome, fallback=True)
                else:
                    # if the recipe is "fasta" we first initialiaze the genome, based on the provided path to the input FASTA file
                    genome, _ = rgc.initialize_genome(
                        fasta_path=specified_args["fasta"],
                        alias=ori_genome,
                        skip_alias_write=True,
                    )
            else:
                try:
                    genome = rgc.get_genome_alias_digest(genome, fallback=True)
                except UndefinedAliasError:
                    _LOGGER.error(
                        "Genome '{}' has not been initialized yet; "
                        "no key found for this alias".format(genome)
                    )
                    return False
            recipe_name = None
            genome_outfolder = os.path.join(rgc.data_dir, genome)
            build_stats_dir = os.path.abspath(
                os.path.join(genome_outfolder, asset_key, asset_tag, BUILD_STATS_DIR)
            )
            locked_map_gencfg = os.path.join(build_stats_dir, LOCKED_BUILD_MAP_CFG)
            map_gencfg = os.path.join(build_stats_dir, BUILD_MAP_CFG)
            is_built, rgc_map = _build_asset(
                genome,
                asset_key,
                asset_tag,
                asset_build_package,
                genome_outfolder,
                specified_args,
                specified_params,
                ori_genome,
                **input_assets,
            )
            if not is_built:
                log_path = os.path.abspath(
                    os.path.join(
                        genome_outfolder,
                        asset_key,
                        asset_tag,
                        BUILD_STATS_DIR,
                        ORI_LOG_NAME,
                    )
                )
                _LOGGER.info(
                    "'{}/{}:{}' was not added to the config, but directory has been left in place. "
                    "See the log file for details: {}".format(
                        genome, asset_key, asset_tag, log_path
                    )
                )
                return False
            _LOGGER.info("Finished building '{}' asset".format(asset_key))
            with rgc_map as r:
                # update asset relationships
                r.update_relatives_assets(
                    genome, asset_key, asset_tag, parent_assets
                )  # adds parents
                for i in parent_assets:
                    parsed_parent = parse_registry_path(i)
                    # adds child (currently built asset) to the parent
                    r.update_relatives_assets(
                        parsed_parent["genome"],
                        parsed_parent["asset"],
                        parsed_parent["tag"],
                        ["{}/{}:{}".format(genome, asset_key, asset_tag)],
                        True,
                    )
                if args.genome_description is not None:
                    _LOGGER.debug(
                        "adding genome ({}) description: '{}'".format(
                            genome, args.genome_description
                        )
                    )
                    r.update_genomes(
                        genome, {CFG_GENOME_DESC_KEY: args.genome_description}
                    )
                if args.tag_description is not None:
                    _LOGGER.debug(
                        "adding tag ({}/{}:{}) description: '{}'".format(
                            genome, asset_key, asset_tag, args.tag_description
                        )
                    )
                    r.update_tags(
                        genome,
                        asset_key,
                        asset_tag,
                        {CFG_TAG_DESC_KEY: args.tag_description},
                    )
            rgc_map._symlink_alias(genome, asset_key, asset_tag)
            if args.map:
                # move the contents of the locked map config to a map config,
                # which is discoverable by the reduce step
                os.rename(locked_map_gencfg, map_gencfg)
                _LOGGER.info(
                    f"Asset metadata saved in '{map_gencfg}'. "
                    f"To make the asset accessible globally run 'refgenie build --reduce'"
                )
        else:
            _raise_missing_recipe_error(recipe_name)

        return True


def _handle_sigint(gat):
    """
    SIGINT handler, unlocks the config file and exists the program

    :param list gat: a list of genome, asset and tag. Used for a message generation.
    :return function: the SIGINT handling function
    """

    def handle(sig, frame):
        _LOGGER.warning("\nThe build was interrupted: {}/{}:{}".format(*gat))
        sys.exit(0)

    return handle


def _check_recipe(recipe):
    """
    Check whether there are any key name clashes in the recipe requirements
    and raise an error if there are

    :param dict recipe: asset_build_package
    :raise ValueError: if any key names are duplicated
    """
    # experimental feature; recipe jsonschema validation
    from jsonschema import validate
    from yacman import load_yaml

    SCHEMA_SRC = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "schemas", "recipe_schema.yaml"
    )
    if os.path.exists(SCHEMA_SRC):
        validate(recipe, load_yaml(filepath=SCHEMA_SRC))
        _LOGGER.info(
            "Recipe validated successfully against a schema: {}".format(SCHEMA_SRC)
        )
    else:
        _LOGGER.warning("Recipe schema not found: {}".format(SCHEMA_SRC))
    # end of validation
    req_keys = []
    for req in [REQ_PARAMS, REQ_ASSETS, REQ_FILES]:
        req_keys.extend([req_dict[KEY] for req_dict in recipe[req]])
    unique = []
    for k in req_keys:
        if k not in unique:
            unique.append(k)
        else:
            raise ValueError(
                "The recipe contains a duplicated requirement"
                " key '{}', which is not permitted.".format(k)
            )
    return recipe


def _seek(
    rgc, genome_name, asset_name, tag_name=None, seek_key=None, enclosing_dir=False
):
    """
    Strict seek. Most use cases in this package require file existence
     check in seek. This function makes it easier
    """
    return rgc.seek_src(
        genome_name=genome_name,
        asset_name=asset_name,
        tag_name=tag_name,
        seek_key=seek_key,
        enclosing_dir=enclosing_dir,
        strict_exists=True,
    )
