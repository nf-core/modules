#!/usr/bin/env python3
"""
Run cellranger aggr with a pre-built aggregation CSV.

Copyright (c) Alena Boos 2026 - MIT License
"""

from subprocess import run

prefix = "${task.ext.prefix ?: meta.id}"
aggr_csv = "${aggr_csv}"
args = "${task.ext.args ?: ''}"

# fmt: off
run(
    [
        "cellranger", "aggr",
        "--id", prefix,
        "--csv", aggr_csv,
        "--localcores", "${task.cpus}",
        "--localmem", "${task.memory.toGiga()}",
    ] + (args.split() if args else []),
    check=True,
)
# fmt: on

version = (
    run(
        ["cellranger", "--version"],
        text=True,
        check=True,
        capture_output=True,
    )
    .stdout.strip()
    .replace("cellranger cellranger-", "")
)

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    cellranger: "{version}"\\n')
