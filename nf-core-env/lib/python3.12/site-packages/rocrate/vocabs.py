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

import sys
import json
if sys.version_info.minor < 9:
    import pkg_resources
else:
    import importlib.resources

# FIXME: Avoid eager loading?
if sys.version_info.minor < 9:
    RO_CRATE = json.loads(pkg_resources.resource_string(
        __name__, "data/ro-crate.jsonld"
    ))
    SCHEMA = json.loads(pkg_resources.resource_string(
        __name__, "data/schema.jsonld"
    ))
else:
    RO_CRATE = json.loads(
        importlib.resources.files(__package__).joinpath("data/ro-crate.jsonld").read_text("utf8")
    )
    SCHEMA = json.loads(
        importlib.resources.files(__package__).joinpath("data/schema.jsonld").read_text("utf8")
    )
SCHEMA_MAP = dict((e["@id"], e) for e in SCHEMA["@graph"])


def term_to_uri(name):
    # NOTE: Assumes RO-Crate's flat-style context
    return RO_CRATE["@context"][name]


def schema_doc(uri):
    # NOTE: Ensure rdfs:comment still appears in newer schema.org downloads
    # TODO: Support terms outside schema.org?
    return SCHEMA_MAP[uri].get("rdfs:comment", "")
