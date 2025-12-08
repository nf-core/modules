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

"""\
Generate a Workflow Testing RO-Crate from a "best practices" Galaxy
workflow repository.

https://github.com/galaxyproject/iwc/blob/main/workflows/README.md
"""

import json
import os
import warnings
from pathlib import Path

import yaml
from rocrate.model.entity import Entity

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
from .common import CrateBuilder


DOCKSTORE_CONF_BASENAME = ".dockstore.yml"
PLANEMO_TEST_SUFFIXES = ["-tests", "_tests", "-test", "_test"]
PLANEMO_TEST_EXTENSIONS = [".yml", ".yaml", ".json"]


def find_workflow(root_dir):
    candidates = [_.name for _ in os.scandir(root_dir) if _.is_file() and _.name.endswith(".ga")]
    if not candidates:
        raise RuntimeError("workflow (.ga file) not found")
    if len(candidates) > 1:
        warnings.warn("Multiple .ga files found, picking one arbitrarily")
        for c in candidates:
            if c.lower() == "main.ga":
                return root_dir / c
    return root_dir / candidates[0]


def find_test_definition(root_dir, wf_name):
    tag = os.path.splitext(wf_name)[0]
    for suffix in PLANEMO_TEST_SUFFIXES:
        for ext in PLANEMO_TEST_EXTENSIONS:
            def_path = root_dir / f"{tag}{suffix}{ext}"
            if def_path.is_file():
                return def_path


def get_workflow_name(root_dir, workflow_relpath):
    dockstore_conf_path = root_dir / DOCKSTORE_CONF_BASENAME
    if not dockstore_conf_path.is_file():
        return None
    try:
        with open(dockstore_conf_path) as f:
            conf = yaml.load(f, Loader=Loader)
    except (OSError, yaml.error.YAMLError):
        return None
    for record in conf.get("workflows", []):
        relpath = record.get("primaryDescriptorPath", "").strip().lstrip("/")
        if Path(relpath) == workflow_relpath:
            return f"{root_dir.name}/{record.get('name')}"


class GalaxyCrateBuilder(CrateBuilder):

    CI_WORKFLOW = "wftest.yml"
    DATA_ENTITIES = [
        ("test-data", "Dataset", "Data files for testing the workflow"),
        ("CHANGELOG.md", "File", "History of changes made to the workflow"),
        ("README.md", "File", "Workflow documentation"),
        (DOCKSTORE_CONF_BASENAME, "File", "Dockstore metadata file"),
        # Planemo test file is treated specially (TestDefinition)
    ]

    @property
    def lang(self):
        return "galaxy"

    def __add_creator(self, workflow, wf_code):
        for c in wf_code.get("creator", []):
            type_ = c.get("class")
            if type_ not in {"Organization", "Person"}:
                continue
            properties = {"@type": type_}
            name = c.get("name")
            if name:
                properties["name"] = name
            id_ = c.get("identifier") if type_ == "Person" else c.get("url")
            creator = self.crate.add(Entity(self.crate, identifier=id_, properties=properties))
            workflow.append_to("creator", creator)

    def add_workflow(
        self,
        wf_source,
        wf_name=None,
        wf_version=None,
        lang_version=None,
        license=None,
        diagram=None,
    ):
        if not wf_name:
            wf_name = get_workflow_name(self.root, wf_source.relative_to(self.root))
        workflow = super().add_workflow(
            wf_source,
            wf_name=wf_name,
            wf_version=wf_version,
            lang_version=lang_version,
            license=license,
            diagram=diagram,
        )
        with open(wf_source) as f:
            try:
                wf_code = json.load(f)
            except json.decoder.JSONDecodeError:
                wf_code = {}
        if "release" in wf_code:
            workflow.setdefault("version", wf_code["release"])
        if "license" in wf_code:
            self.crate.root_dataset.setdefault("license", wf_code["license"])
        self.__add_creator(workflow, wf_code)
        return workflow

    def add_test_suite(self, workflow=None, ci_workflow=None):
        suite = super().add_test_suite(workflow, ci_workflow=ci_workflow)
        if workflow is None:
            workflow = suite["mainEntity"]
        def_path = find_test_definition(self.root, workflow.id)
        if def_path:
            self.crate.add_test_definition(
                suite, source=def_path, dest_path=def_path.relative_to(self.root), engine="planemo"
            )
        return suite


def make_crate(
    root,
    workflow=None,
    repo_url=None,
    wf_name=None,
    wf_version=None,
    lang_version=None,
    license=None,
    ci_workflow=None,
    diagram=None,
):
    builder = GalaxyCrateBuilder(root, repo_url=repo_url)
    if not workflow:
        workflow = find_workflow(root)
    return builder.build(
        workflow,
        wf_name=wf_name,
        wf_version=wf_version,
        lang_version=lang_version,
        license=license,
        ci_workflow=ci_workflow,
        diagram=diagram,
    )
