# Copyright 2022 CRS4.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.

from pathlib import Path

import click
from . import find_workflow, LANG_MODULES, __version__


@click.command()
@click.option(
    "-r",
    "--root",
    type=click.Path(exists=True, file_okay=False, readable=True, path_type=Path),
    help="workflow repository root",
    default=Path.cwd,
)
@click.option(
    "-l",
    "--lang",
    type=click.Choice(list(LANG_MODULES)),
    help="workflow language (default: auto-detect)",
)
@click.option(
    "-w", "--workflow", type=click.Path(path_type=Path), help="workflow file (default: auto-detect)"
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help=(
        "output directory or zip file. The default is the repository root itself, in which case"
        " only the metadata file is written"
    ),
)
@click.option("--repo-url", help="workflow repository URL")
@click.option("--wf-name", help="workflow name")
@click.option("--wf-version", help="workflow version")
@click.option("--lang-version", help="workflow language version")
@click.option("--license", help="license URL")
@click.option(
    "--ci-workflow",
    help="filename (basename) of the GitHub Actions workflow that runs the tests for the workflow",
)
@click.option("--diagram", help="relative path of the workflow diagram")
@click.option("--version", help="print version and exit", is_flag=True)
def cli(
    root,
    lang,
    workflow,
    output,
    repo_url,
    wf_name,
    wf_version,
    lang_version,
    license,
    ci_workflow,
    diagram,
    version,
):
    if version:
        print(__version__)
        return
    auto_workflow = None
    if not lang:
        lang, auto_workflow = find_workflow(root)
    if not workflow:
        workflow = auto_workflow
    make_crate = LANG_MODULES[lang].make_crate
    crate = make_crate(
        root,
        workflow=workflow,
        repo_url=repo_url,
        wf_name=wf_name,
        wf_version=wf_version,
        lang_version=lang_version,
        license=license,
        ci_workflow=ci_workflow,
        diagram=diagram,
    )
    if not output:
        output = root
    try:
        in_place = output.samefile(root)
    except OSError:
        in_place = False
    if in_place:
        crate.metadata.write(output)
    elif output.suffix == ".zip":
        crate.write_zip(output)
    else:
        crate.write(output)


if __name__ == "__main__":
    cli()
