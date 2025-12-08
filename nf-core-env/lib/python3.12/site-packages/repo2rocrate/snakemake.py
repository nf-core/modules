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
Generate a Workflow Testing RO-Crate from a "best practices" Snakemake
workflow repository.

https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
https://snakemake.github.io/snakemake-workflow-catalog/?rules=true
"""

import re

from .common import CrateBuilder


WF_BASENAME = "Snakefile"


def find_workflow(root_dir):
    candidates = [
        root_dir / "workflow" / WF_BASENAME,
        root_dir / WF_BASENAME,
    ]
    for p in candidates:
        if p.is_file():
            return p
    raise RuntimeError(f"workflow definition (one of: {', '.join(map(str, candidates))}) not found")


def get_lang_version(workflow_path):
    pattern = re.compile(r"""min_version\(\s*['"]([^'"]+)['"]\s*\)""")
    with open(workflow_path) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                return m.groups()[0]


class SnakemakeCrateBuilder(CrateBuilder):

    DATA_ENTITIES = [
        (".tests/integration", "Dataset", "Integration tests for the workflow"),
        (".tests/unit", "Dataset", "Unit tests for the workflow"),
        ("workflow", "Dataset", "Workflow folder"),
        ("config", "Dataset", "Configuration folder"),
        ("results", "Dataset", "Output files"),
        ("resources", "Dataset", "Retrieved resources"),
        ("workflow/rules", "Dataset", "Workflow rule modules"),
        ("workflow/envs", "Dataset", "Conda environments"),
        ("workflow/scripts", "Dataset", "Scripts folder"),
        ("workflow/notebooks", "Dataset", "Notebooks folder"),
        ("workflow/report", "Dataset", "Report caption files"),
        ("workflow/schemas", "Dataset", "Validation files"),
    ]
    for s in "README", "LICENSE", "CODE_OF_CONDUCT", "CONTRIBUTING":
        for ext in "", ".md":
            DATA_ENTITIES.append((s + ext, "File", ""))
    DIAGRAM = "images/rulegraph.svg"

    @property
    def lang(self):
        return "snakemake"

    def add_workflow(
        self,
        wf_source,
        wf_name=None,
        wf_version=None,
        lang_version=None,
        license=None,
        diagram=None,
    ):
        if not lang_version:
            lang_version = get_lang_version(wf_source)
        workflow = super().add_workflow(
            wf_source,
            wf_name=wf_name,
            wf_version=wf_version,
            lang_version=lang_version,
            license=license,
            diagram=diagram,
        )
        return workflow


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
    builder = SnakemakeCrateBuilder(root, repo_url=repo_url)
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
