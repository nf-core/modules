# Copyright 2019-2025 The University of Manchester, UK
# Copyright 2020-2025 Vlaams Instituut voor Biotechnologie (VIB), BE
# Copyright 2020-2025 Barcelona Supercomputing Center (BSC), ES
# Copyright 2020-2025 Center for Advanced Studies, Research and Development in Sardinia (CRS4), IT
# Copyright 2022-2025 École Polytechnique Fédérale de Lausanne, CH
# Copyright 2024-2025 Data Centre, SciLifeLab, SE
# Copyright 2024-2025 National Institute of Informatics (NII), JP
# Copyright 2025 Senckenberg Society for Nature Research (SGN), DE
# Copyright 2025 European Molecular Biology Laboratory (EMBL), Heidelberg, DE
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from pathlib import Path

import click
from .rocrate import ROCrate
from .model.computerlanguage import LANG_MAP
from .model.testservice import SERVICE_MAP
from .model.softwareapplication import APP_MAP


LANG_CHOICES = list(LANG_MAP)
SERVICE_CHOICES = list(SERVICE_MAP)
ENGINE_CHOICES = list(APP_MAP)


class CSVParamType(click.ParamType):
    name = "csv"

    def convert(self, value, param, ctx):
        if isinstance(value, (list, tuple, set, frozenset)):
            return value
        try:
            return value.split(",") if value else []
        except AttributeError:
            self.fail(f"{value!r} is not splittable", param, ctx)


class KeyValueParamType(click.ParamType):
    name = "key_value"

    def convert(self, value, param, ctx):
        try:
            return tuple(value.split("=", 1)) if value else ()
        except AttributeError:
            self.fail(f"{value!r} is not splittable", param, ctx)


CSV = CSVParamType()
KeyValue = KeyValueParamType()
OPTION_CRATE_PATH = click.option(
    "-c",
    "--crate-dir",
    type=click.Path(),
    default=os.getcwd,
    help="The path to the root data entity of the crate. Defaults to the current working directory.",
)
OPTION_PROPS = click.option(
    "-P",
    "--property",
    type=KeyValue,
    multiple=True,
    metavar="KEY=VALUE",
    help="Add an additional property to the metadata for this entity. Can be used multiple times to set multiple properties.",
)


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--gen-preview", is_flag=True, help="Generate a HTML preview file for the crate."
)
@click.option(
    "-e",
    "--exclude",
    type=CSV,
    metavar="NAME",
    help="Exclude files or directories from the metadata file. NAME may be a single name or a comma-separated list of names.",
)
@OPTION_CRATE_PATH
def init(crate_dir, gen_preview, exclude):
    crate = ROCrate(crate_dir, init=True, gen_preview=gen_preview, exclude=exclude)
    crate.metadata.write(crate_dir)
    if crate.preview:
        crate.preview.write(crate_dir)


@cli.group()
def add():
    pass


@add.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@OPTION_CRATE_PATH
@OPTION_PROPS
def file(crate_dir, path, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    source = Path(path).resolve(strict=True)
    try:
        dest_path = source.relative_to(crate_dir)
    except ValueError:
        # For now, only support adding an existing file to the metadata
        raise ValueError(f"{source} is not in the crate dir {crate_dir}")
    crate.add_file(source, dest_path, properties=dict(property))
    crate.metadata.write(crate_dir)


@add.command()
@click.argument("path", type=click.Path(exists=True, file_okay=False))
@OPTION_CRATE_PATH
@OPTION_PROPS
def dataset(crate_dir, path, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    source = Path(path).resolve(strict=True)
    try:
        dest_path = source.relative_to(crate_dir)
    except ValueError:
        # For now, only support adding an existing directory to the metadata
        raise ValueError(f"{source} is not in the crate dir {crate_dir}")
    crate.add_dataset(source, dest_path, properties=dict(property))
    crate.metadata.write(crate_dir)


@add.command()
@click.argument("path", type=click.Path(exists=True))
@click.option(
    "-l",
    "--language",
    type=click.Choice(LANG_CHOICES),
    default="cwl",
    help="The workflow language.",
)
@OPTION_CRATE_PATH
@OPTION_PROPS
def workflow(crate_dir, path, language, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    source = Path(path).resolve(strict=True)
    try:
        dest_path = source.relative_to(crate_dir)
    except ValueError:
        # For now, only support marking an existing file as a workflow
        raise ValueError(f"{source} is not in the crate dir {crate_dir}")
    # TODO: add command options for main and gen_cwl
    crate.add_workflow(
        source,
        dest_path,
        main=True,
        lang=language,
        gen_cwl=False,
        properties=dict(property),
    )
    crate.metadata.write(crate_dir)


@add.command(name="test-suite")
@click.option("-i", "--identifier")
@click.option("-n", "--name")
@click.option("-m", "--main-entity")
@OPTION_CRATE_PATH
@OPTION_PROPS
def suite(crate_dir, identifier, name, main_entity, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    suite = crate.add_test_suite(
        identifier=identifier,
        name=name,
        main_entity=main_entity,
        properties=dict(property),
    )
    crate.metadata.write(crate_dir)
    print(suite.id)


@add.command(name="test-instance")
@click.argument("suite")
@click.argument("url")
@click.option("-r", "--resource", default="")
@click.option("-s", "--service", type=click.Choice(SERVICE_CHOICES), default="jenkins")
@click.option("-i", "--identifier")
@click.option("-n", "--name")
@OPTION_CRATE_PATH
@OPTION_PROPS
def instance(crate_dir, suite, url, resource, service, identifier, name, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    instance_ = crate.add_test_instance(
        suite,
        url,
        resource=resource,
        service=service,
        identifier=identifier,
        name=name,
        properties=dict(property),
    )
    crate.metadata.write(crate_dir)
    print(instance_.id)


@add.command(name="test-definition")
@click.argument("suite")
@click.argument("path", type=click.Path(exists=True))
@click.option("-e", "--engine", type=click.Choice(ENGINE_CHOICES), default="planemo")
@click.option("-v", "--engine-version")
@OPTION_CRATE_PATH
@OPTION_PROPS
def definition(crate_dir, suite, path, engine, engine_version, property):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    source = Path(path).resolve(strict=True)
    try:
        dest_path = source.relative_to(crate_dir)
    except ValueError:
        # For now, only support marking an existing file as a test definition
        raise ValueError(f"{source} is not in the crate dir {crate_dir}")
    crate.add_test_definition(
        suite,
        source=source,
        dest_path=dest_path,
        engine=engine,
        engine_version=engine_version,
        properties=dict(property),
    )
    crate.metadata.write(crate_dir)


@cli.command()
@click.argument("dst", type=click.Path(writable=True))
@OPTION_CRATE_PATH
def write_zip(crate_dir, dst):
    crate = ROCrate(crate_dir, init=False, gen_preview=False)
    crate.write_zip(dst)


if __name__ == "__main__":
    cli()
