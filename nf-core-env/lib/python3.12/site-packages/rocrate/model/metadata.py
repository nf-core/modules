#!/usr/bin/env python

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
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
from pathlib import Path

from .file import File
from .dataset import Dataset


WORKFLOW_PROFILE = "https://w3id.org/workflowhub/workflow-ro-crate/1.0"


class Metadata(File):
    """\
    RO-Crate metadata file.
    """
    BASENAME = "ro-crate-metadata.json"
    PROFILE = "https://w3id.org/ro/crate/1.1"

    def __init__(self, crate, source=None, dest_path=None, properties=None):
        if source is None and dest_path is None:
            dest_path = self.BASENAME
        super().__init__(
            crate,
            source=source,
            dest_path=dest_path,
            fetch_remote=False,
            validate_url=False,
            properties=properties
        )
        # https://www.researchobject.org/ro-crate/1.1/appendix/jsonld.html#extending-ro-crate
        self.extra_contexts = []
        self.extra_terms = {}

    def _empty(self):
        # default properties of the metadata entry
        val = {"@id": self.id,
               "@type": "CreativeWork",
               "conformsTo": {"@id": self.PROFILE},
               "about": {"@id": "./"}}
        return val

    # Generate the crate's `ro-crate-metadata.json`.
    # @return [String] The rendered JSON-LD as a "prettified" string.
    def generate(self):
        graph = []
        for entity in self.crate.get_entities():
            graph.append(entity.properties())
        context = [f'{self.PROFILE}/context']
        context.extend(self.extra_contexts)
        if self.extra_terms:
            context.append(self.extra_terms)
        if len(context) == 1:
            context = context[0]
        return {'@context': context, '@graph': graph}

    def stream(self, chunk_size=8192):
        content = self.generate()
        yield self.id, str.encode(json.dumps(content, indent=4, sort_keys=True), encoding='utf-8')

    def _has_writeable_stream(self):
        return True

    def write(self, dest_base):
        write_path = Path(dest_base) / self.id
        super()._write_from_stream(write_path)

    @property
    def root(self) -> Dataset:
        return self.crate.root_dataset


class LegacyMetadata(Metadata):

    BASENAME = "ro-crate-metadata.jsonld"
    PROFILE = "https://w3id.org/ro/crate/1.0"


# https://github.com/ResearchObject/ro-terms/tree/master/test
TESTING_EXTRA_TERMS = {
    "TestSuite": "https://w3id.org/ro/terms/test#TestSuite",
    "TestInstance": "https://w3id.org/ro/terms/test#TestInstance",
    "TestService": "https://w3id.org/ro/terms/test#TestService",
    "TestDefinition": "https://w3id.org/ro/terms/test#TestDefinition",
    "PlanemoEngine": "https://w3id.org/ro/terms/test#PlanemoEngine",
    "JenkinsService": "https://w3id.org/ro/terms/test#JenkinsService",
    "TravisService": "https://w3id.org/ro/terms/test#TravisService",
    "GithubService": "https://w3id.org/ro/terms/test#GithubService",
    "instance": "https://w3id.org/ro/terms/test#instance",
    "runsOn": "https://w3id.org/ro/terms/test#runsOn",
    "resource": "https://w3id.org/ro/terms/test#resource",
    "definition": "https://w3id.org/ro/terms/test#definition",
    "engineVersion": "https://w3id.org/ro/terms/test#engineVersion"
}


def metadata_class(descriptor_id):
    basename = descriptor_id.rsplit("/", 1)[-1]
    if basename == Metadata.BASENAME:
        return Metadata
    elif basename == LegacyMetadata.BASENAME:
        return LegacyMetadata
    else:
        raise ValueError(f"Invalid metadata descriptor ID: {descriptor_id!r}")
