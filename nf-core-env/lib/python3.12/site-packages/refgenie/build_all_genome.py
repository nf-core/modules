"""
A helper script to create SLURM submission scripts for all the assets
defined in asset_build_packages for a given genome
"""
import argparse
import os

import divvy
from ubiquerg import expandpath

from .asset_build_packages import asset_build_packages

parser = argparse.ArgumentParser(
    description="Builds submission scripts for all assets for a genome"
)
parser.add_argument(
    "-g",
    "--genome",
    dest="genome",
    type=str,
    help="genome to build the submission scripts for",
)
parser.add_argument(
    "-p",
    "--path",
    dest="path",
    type=str,
    help="path to the desired submission directory location",
)
parser.add_argument(
    "-pt",
    "--partition",
    dest="PARTITION",
    type=str,
    help="partition in SLURM submission script",
    default="standard",
)
parser.add_argument(
    "-m",
    "--mem",
    dest="MEM",
    type=str,
    help="mem in SLURM submission script",
    default="200000",
)
parser.add_argument(
    "-t",
    "--time",
    dest="TIME",
    type=str,
    help="time in SLURM submission script",
    default="10:00:00",
)
parser.add_argument(
    "-c",
    "--cores",
    dest="CORES",
    type=str,
    help="cpus-per-task in SLURM submission script",
    default="4",
)
parser.add_argument(
    "-o",
    "--output",
    dest="LOGFILE",
    type=str,
    help="output in SLURM submission script",
    default=None,
)
parser.add_argument(
    "-j",
    "--job-name",
    dest="JOBNAME",
    type=str,
    help="job-name in SLURM submission script",
    default=None,
)

args = parser.parse_args()


def _make_sub_dir(path, genome):
    """
    create submission scripts directory

    :param str path: path where the submission scripts directory should be created
    :param str genome: name of the genome
    :return str: created path
    """
    path = os.path.join(expandpath(path), "submission_scripts", genome)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def _req_input_to_args(req_input):
    """
    Given a list of the required inputs for the build command, create an args
    string

    :param list[str] req_input: input names
    :return str: args string
    """
    return ["--" + x + " <arg_here>" for x in req_input]


subdir_path = _make_sub_dir(args.path, args.genome)
dcc = divvy.ComputingConfiguration()
dcc.activate_package("slurm")
cmd_template = "refgenie build -g {g} -a {a} {req_input_str}"
genome = args.genome
to_remove = ["genome", "path"]

data = vars(args)
for i in to_remove:
    data.pop(i)

for asset in asset_build_packages:
    sub_script = os.path.join(subdir_path, asset + ".sub")
    req_input = asset_build_packages[asset]["required_inputs"]
    if req_input:
        print(
            "{} asset requires additional input in the command ({}), so '{}'"
            " requires manual edit".format(asset, req_input, sub_script)
        )
        req_str = " ".join(_req_input_to_args(req_input))
    else:
        req_str = ""
    data["CODE"] = cmd_template.format(g=genome, a=asset, req_input_str=req_str)
    data["LOGFILE"] = asset + ".log"
    data["JOBNAME"] = asset + "Build"
    dcc.write_script(sub_script, data)
