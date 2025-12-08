#!/usr/bin/env python
"""
Each iGenome has the following nested directory structure:
    Species/
    Source/
    Build/
    Annotation/ Sequence/
"""
import argparse
import os
import sys
import tarfile
from glob import glob
from shutil import move

import refgenconf
from refgenconf import get_dir_digest
from refgenconf.const import *
from ubiquerg import mkabs, query_yes_no, untar

from .exceptions import MissingGenomeConfigError
from .refgenie import _seek


def build_argparser():
    """
    Build a parser for this tool

    :return argparse.ArgumentParser: constructed parser
    """
    parser = argparse.ArgumentParser(
        description="Integrates every asset from the downloaded iGenomes"
        " tarball/directory with Refgenie asset management system"
    )
    parser.add_argument(
        "-p",
        "--path",
        dest="path",
        type=str,
        help="path to the desired genome tarball or directory to integrate",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        dest="genome",
        type=str,
        help="name to be assigned to the selected genome",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--config",
        dest="config",
        type=str,
        help="path to local genome configuration file. Optional if '{}' environment variable is set.".format(
            ", ".join(refgenconf.CFG_ENV_VARS)
        ),
        required=False,
    )
    return parser


def untar_or_copy(p, dest):
    """
    Depending on a kind of the provided path, either copy or extract it to the destination directory

    :param str p: path to the directory to be copied or tarball to be extracted
    :param str dest: where to extract file or copy dir
    :return bool: whether the process was successful
    """
    if os.path.exists(p):
        if os.path.isdir(p):
            fun = move
            dest = os.path.join(dest, os.path.basename(p))
        elif tarfile.is_tarfile(p):
            fun = untar
            print("Extracting '{}'".format(p))
        else:
            raise ValueError("Provided path is neither a directory nor a tar archive.")
        fun(p, dest)
        print("Moved '{}' to '{}'".format(p, dest))
        return True
    return False


def refgenie_add(rgc, asset_dict, path, force=False):
    """
    Add an external asset to the config.
    File existence is checked and asset files are transferred to the selected
    tag subdirectory

    :param refgenconf.RefGenConf rgc: genome configuration object
    :param dict asset_dict: a single parsed registry path
    :param str path: the path provided by the user. Must be relative to the
        specific genome directory
    :param bool force: whether the replacement of a possibly existing asset
        should be forced
    """
    # remove the first directory from the provided path if it is the genome name
    path = (
        os.path.join(*path.split(os.sep)[1:])
        if path.split(os.sep)[0] == asset_dict["genome"]
        else path
    )
    tag = asset_dict["tag"] or rgc.get_default_tag(
        asset_dict["genome"], asset_dict["asset"]
    )
    outfolder = os.path.abspath(os.path.join(rgc[CFG_FOLDER_KEY], asset_dict["genome"]))
    abs_asset_path = os.path.join(outfolder, path)
    if asset_dict["seek_key"] is None:
        # if seek_key is not specified we're about to move a directory to
        # the tag subdir
        tag_path = os.path.join(abs_asset_path, tag)
        from shutil import copytree as cp
    else:
        # if seek_key is specified we're about to move just a single file to
        # he tag subdir
        tag_path = os.path.join(os.path.dirname(abs_asset_path), tag)
        if not os.path.exists(tag_path):
            os.makedirs(tag_path)
        from shutil import copy2 as cp
    if os.path.exists(abs_asset_path):
        if not os.path.exists(tag_path):
            cp(abs_asset_path, tag_path)
        else:
            if not force and not query_yes_no(
                "Path '{}' exists. Do you want to overwrite?".format(tag_path)
            ):
                return False
            else:
                _remove(tag_path)
                cp(abs_asset_path, tag_path)
    else:
        raise OSError(
            "Absolute path '{}' does not exist. "
            "The provided path must be relative to: {}".format(
                abs_asset_path, rgc[CFG_FOLDER_KEY]
            )
        )
    rgc.make_writable()
    gat_bundle = [asset_dict["genome"], asset_dict["asset"], tag]
    td = {
        CFG_ASSET_PATH_KEY: path
        if os.path.isdir(abs_asset_path)
        else os.path.dirname(path)
    }
    rgc.update_tags(*gat_bundle, data=td)
    # seek_key points to the entire dir if not specified
    seek_key_value = (
        os.path.basename(abs_asset_path) if asset_dict["seek_key"] is not None else "."
    )
    sk = {asset_dict["seek_key"] or asset_dict["asset"]: seek_key_value}
    rgc.update_seek_keys(*gat_bundle, keys=sk)
    rgc.set_default_pointer(asset_dict["genome"], asset_dict["asset"], tag)
    # a separate update_tags call since we want to use the get_asset method
    # that requires a complete asset entry in rgc
    td = {CFG_ASSET_CHECKSUM_KEY: get_dir_digest(_seek(rgc, *gat_bundle))}
    rgc.update_tags(*gat_bundle, data=td)
    # Write the updated refgenie genome configuration
    rgc.write()
    rgc.make_readonly()
    return True


def main():
    """main workflow"""
    parser = build_argparser()
    args, remaining_args = parser.parse_known_args()
    cfg = refgenconf.select_genome_config(
        filename=args.config, check_exist=True, strict_env=True
    )
    if not cfg:
        raise MissingGenomeConfigError(args.config)
    rgc = refgenconf.RefGenConf(filepath=cfg, writable=False)
    pths = [args.path, mkabs(args.path, rgc.genome_folder)]
    if not untar_or_copy(
        pths[0], os.path.join(rgc.genome_folder, args.genome)
    ) and not untar_or_copy(pths[1], os.path.join(rgc.genome_folder, args.genome)):
        raise OSError(
            "Path '{}' does not exist. Tried: {}".format(args.path, " and ".join(pths))
        )
    path_components = [rgc.genome_folder] + [args.genome] + ["*"] * 3 + ["Sequence"]
    assets_paths = glob(os.path.join(*path_components))
    assert len(assets_paths) > 0, OSError(
        "Your iGenomes directory is corrupted, more than one directory matched by {}."
        "\nMatched dirs: {}".format(
            os.path.join(*path_components), ", ".join(assets_paths)
        )
    )
    assets_path = assets_paths[0]
    asset_names = [d for d in os.listdir(assets_path) if os.path.isdir(assets_path)]
    processed = []
    for a in asset_names:
        asset_dict = {"genome": args.genome, "asset": a, "tag": None, "seek_key": None}
        asset_path = os.path.relpath(os.path.join(assets_path, a), rgc.genome_folder)
        if refgenie_add(rgc, asset_dict, asset_path):
            processed.append("{}/{}".format(asset_dict["genome"], asset_dict["asset"]))
    print("Added assets: \n- {}".format("\n- ".join(processed)))


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
        raise ValueError("path '{}' is neither a file nor a dir.".format(path))
    return path


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
