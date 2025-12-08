"""
Update a nextflow.config file with refgenie genomes
"""

import logging
import os
import re
from pathlib import Path
from textwrap import dedent

import rich
import rich.traceback
import yaml

import nf_core.utils

# import refgenconf


# Set up logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

# # Setup rich traceback
stderr = rich.console.Console(stderr=True, force_terminal=nf_core.utils.rich_force_colors())
rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)

NF_CFG_TEMPLATE = """
// This is a read-only config file managed by refgenie. Manual changes to this file will be overwritten
// To make changes here, use refgenie to update the reference genome data
params {{
  genomes {{
{content}
  }}
}}
"""


def _print_nf_config(rgc):
    """
    Generate a nextflow config file with the genomes
    from the refgenie config file
    Adapted from: https://github.com/refgenie/refgenie_nfcore

    Takes a RefGenConf object as argument
    """
    abg = rgc.list_assets_by_genome()
    genomes_str = ""
    alias_translations = _get_alias_translation_file(rgc)
    for genome, asset_list in abg.items():
        genomes_str += f"    '{genome}' {{\n"
        for asset in asset_list:
            try:
                pth = rgc.seek(genome, asset)
            # Catch general exception instead of refgencof exception --> no refgenconf import needed
            except Exception:
                log.warning(f"{genome}/{asset} is incomplete, ignoring...")
            else:
                # Translate an alias name to the alias used in the pipeline
                if asset in alias_translations.keys():
                    log.info(f"Translating refgenie asset alias {asset} to {alias_translations[asset]}.")
                    asset = alias_translations[asset]
                genomes_str += f'      {asset.ljust(20, " ")} = "{pth}"\n'
        genomes_str += "    }\n"

    return NF_CFG_TEMPLATE.format(content=genomes_str)


def _update_nextflow_home_config(refgenie_genomes_config_file, nxf_home):
    """
    Update the $NXF_HOME/config file by adding a includeConfig statement to it
    for the 'refgenie_genomes_config_file' if not already defined
    """
    # Check if NXF_HOME/config exists and has a
    include_config_string = dedent(
        f"""
        ///// >>> nf-core + RefGenie >>> /////
        // !! Contents within this block are managed by 'nf-core/tools' !!
        // Includes auto-generated config file with RefGenie genome assets
        includeConfig '{os.path.abspath(refgenie_genomes_config_file)}'
        ///// <<< nf-core + RefGenie <<< /////
        """
    )
    nxf_home_config = Path(nxf_home) / "config"
    if os.path.exists(nxf_home_config):
        # look for include statement in config
        has_include_statement = False
        with open(nxf_home_config) as fh:
            lines = fh.readlines()
            for line in lines:
                if re.match(rf"\s*includeConfig\s*'{os.path.abspath(refgenie_genomes_config_file)}'", line):
                    has_include_statement = True
                    break

        # if include statement is missing, add it to the last line
        if not has_include_statement:
            with open(nxf_home_config, "a") as fh:
                fh.write(include_config_string)

            log.info(f"Included refgenie_genomes.config to {nxf_home_config}")

    else:
        # create new config and add include statement
        with open(nxf_home_config, "w") as fh:
            fh.write(include_config_string)
            log.info(f"Created new nextflow config file: {nxf_home_config}")


def _get_alias_translation_file(rgc):
    """
    Read a file containing alias translations.

    Alias translation file should be located in the same folder as the refgenie `genome_config.yaml` file,
    the path is set to $REFGENIE environment variable by `refgenie init`.
    Alias translation file should be named `alias_translations.yaml`

    Input file contains the name of refgenie server aliases as keys and the name of the respective nf-core pipeline aliases as values.
    Such as:
    ensembl_gtf: gtf
    star_index: star
    """
    translations = {}

    if "REFGENIE" in os.environ:
        refgenie_genomes_config_path = os.environ.get("REFGENIE")
        refgenie_genomes_config_directory = Path(refgenie_genomes_config_path).parents[0]
    elif "genome_folder" in rgc:
        refgenie_genomes_config_directory = Path(rgc["genome_folder"])
    else:
        return translations

    try:
        with open(refgenie_genomes_config_directory / "alias_translations.yaml") as yaml_file:
            translations = yaml.load(yaml_file, Loader=yaml.Loader)
    except FileNotFoundError:
        pass

    return translations


def update_config(rgc):
    """
    Update the genomes.config file after a local refgenie database has been updated

    This function is executed after running 'refgenie pull <genome>/<asset>'
    The refgenie config file is transformed into a nextflow.config file, which is used to
    overwrite the 'refgenie_genomes.config' file.
    The path to the target  config file is inferred from the following options, in order:

    - the 'nextflow_config' attribute in the refgenie config file
    - the NXF_REFGENIE_PATH environment variable
    - otherwise defaults to: $NXF_HOME/nf-core/refgenie_genomes.config

    Additionally, an 'includeConfig' statement is added to the file $NXF_HOME/config
    """

    # Compile nextflow refgenie_genomes.config from refgenie config
    refgenie_genomes = _print_nf_config(rgc)

    # Get the path to NXF_HOME
    # If NXF_HOME is not set, create it at $HOME/.nextflow
    # If $HOME is not set, set nxf_home to false
    nxf_home = os.environ.get("NXF_HOME")
    if not nxf_home:
        try:
            nxf_home = Path.home() / ".nextflow"
            if not os.path.exists(nxf_home):
                log.info(f"Creating NXF_HOME directory at {nxf_home}")
                os.makedirs(nxf_home, exist_ok=True)
        except RuntimeError:
            nxf_home = False

    # Get the path for storing the updated refgenie_genomes.config
    if hasattr(rgc, "nextflow_config"):
        refgenie_genomes_config_file = rgc.nextflow_config
    elif "NXF_REFGENIE_PATH" in os.environ:
        refgenie_genomes_config_file = os.environ.get("NXF_REFGENIE_PATH")
    elif nxf_home:
        refgenie_genomes_config_file = Path(nxf_home) / "nf-core/refgenie_genomes.config"
    else:
        log.info("Could not determine path to 'refgenie_genomes.config' file.")
        return False

    if refgenie_genomes_config_file is None:
        log.info("Could not determine path to 'refgenie_genomes.config' file.")
        return False

    # Save the updated genome config
    try:
        with open(refgenie_genomes_config_file, "w+") as fh:
            fh.write(refgenie_genomes)
        log.info(f"Updated nf-core genomes config: {refgenie_genomes_config_file}")
    except FileNotFoundError:
        log.info(f"Could not write to {refgenie_genomes_config_file}")
        return False

    # Add include statement to NXF_HOME/config
    if nxf_home:
        _update_nextflow_home_config(refgenie_genomes_config_file, nxf_home)

    return True
