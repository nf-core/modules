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

from abc import ABCMeta, abstractmethod
from pathlib import Path

from rocrate.rocrate import ROCrate
from .utils import get_ci_wf_endpoint

GH_API_URL = "https://api.github.com"


class CrateBuilder(metaclass=ABCMeta):

    CI_WORKFLOW = "main.yml"
    DATA_ENTITIES = []
    DIAGRAM = None

    def __init__(self, root, repo_url=None):
        self.root = Path(root)
        self.repo_url = repo_url
        self.crate = ROCrate(gen_preview=False)

    @property
    @abstractmethod
    def lang(self):
        pass

    def build(
        self,
        wf_source,
        wf_name=None,
        wf_version=None,
        lang_version=None,
        license=None,
        ci_workflow=None,
        diagram=None,
    ):
        workflow = self.add_workflow(
            wf_source,
            wf_name=wf_name,
            wf_version=wf_version,
            lang_version=lang_version,
            license=license,
            diagram=diagram,
        )
        self.add_test_suite(workflow=workflow, ci_workflow=ci_workflow)
        self.add_data_entities()
        return self.crate

    def add_workflow(
        self,
        wf_source,
        wf_name=None,
        wf_version=None,
        lang_version=None,
        license=None,
        diagram=None,
    ):
        if not diagram:
            diagram = self.DIAGRAM
        workflow = self.crate.add_workflow(
            wf_source,
            wf_source.relative_to(self.root),
            main=True,
            lang=self.lang,
            lang_version=lang_version,
            gen_cwl=False,
        )
        workflow["name"] = self.crate.root_dataset["name"] = wf_name or self.root.name
        if wf_version:
            workflow["version"] = wf_version
        # Is there a Python equivalent of https://github.com/licensee/licensee?
        if license:
            self.crate.root_dataset["license"] = license
        if self.repo_url:
            workflow["url"] = self.crate.root_dataset["isBasedOn"] = self.repo_url
        if diagram:
            diag_source = self.root / diagram
            if diag_source.is_file():
                diag = self.crate.add_file(
                    diag_source,
                    diagram,
                    properties={
                        "name": "Workflow diagram",
                        "@type": ["File", "ImageObject"],
                    },
                )
                workflow["image"] = diag
        return workflow

    def add_test_suite(self, workflow=None, ci_workflow=None):
        if not self.repo_url:
            return None
        if not ci_workflow:
            ci_workflow = self.CI_WORKFLOW
        ci_wf_source = self.root / ".github" / "workflows" / ci_workflow
        if not ci_wf_source.is_file():
            return None
        wf_name = workflow["name"] if workflow else self.crate.mainEntity["name"]
        suite = self.crate.add_test_suite(name=f"Test suite for {wf_name}", main_entity=workflow)
        resource = get_ci_wf_endpoint(self.repo_url, ci_workflow)
        self.crate.add_test_instance(
            suite,
            GH_API_URL,
            resource=resource,
            service="github",
            name=f"GitHub Actions workflow for testing {wf_name}",
        )
        return suite

    def add_data_entities(self, data_entities=None):
        if not data_entities:
            data_entities = self.DATA_ENTITIES
        for relpath, type_, description in data_entities:
            if self.crate.get(str(relpath)):
                continue  # existing entity is more specific
            source = self.root / relpath
            properties = {"description": description} if description else None
            if type_ == "File":
                if source.is_file():
                    self.crate.add_file(source, relpath, properties=properties)
            elif type_ == "Dataset":
                if source.is_dir():
                    self.crate.add_dataset(source, relpath, properties=properties)
            else:
                raise ValueError(f"Unexpected type: {type_!r}")
