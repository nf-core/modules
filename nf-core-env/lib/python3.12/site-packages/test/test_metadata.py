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

import pytest
from copy import deepcopy

from rocrate.metadata import find_root_entity_id


@pytest.mark.parametrize("root,basename", [
    ("", "ro-crate-metadata.json"),
    ("", "ro-crate-metadata.jsonld"),
    ("https://example.org/crate/", "ro-crate-metadata.json"),
    ("https://example.org/crate/", "ro-crate-metadata.jsonld"),
    ("", "bad-name.json"),
])
def test_find_root(root, basename):
    metadata_id = root + basename
    root_id = root or "./"
    entities = {_["@id"]: _ for _ in [
        {
            "@id": metadata_id,
            "@type": "CreativeWork",
            "about": {"@id": root_id},
            "conformsTo": [
                {"@id": "https://w3id.org/ro/crate/1.1"},
                {"@id": "https://example.org/fancy-ro-crate/1.0"},
            ]
        },
        {
            "@id": root_id,
            "@type": "Dataset",
        },
    ]}
    if basename not in {"ro-crate-metadata.json", "ro-crate-metadata.jsonld"}:
        with pytest.raises(KeyError):
            find_root_entity_id(entities)
    else:
        assert find_root_entity_id(entities) == (metadata_id, root_id)


def test_find_root_bad_entities():
    orig_entities = {
        "ro-crate-metadata.json": {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
        },
        "./": {
            "@id": "./",
            "@type": "Dataset",
        },
    }
    # missing "about"
    entities = deepcopy(orig_entities)
    del entities["ro-crate-metadata.json"]["about"]
    with pytest.raises(ValueError, match="does not reference"):
        find_root_entity_id(entities)
    # "about" does not reference the root entity
    entities = deepcopy(orig_entities)
    for about in "http://example.org", {"@id": "http://example.org"}:
        entities["ro-crate-metadata.json"]["about"] = about
        with pytest.raises(ValueError, match="does not reference"):
            find_root_entity_id(entities)
    # metadata type is not CreativeWork
    entities = deepcopy(orig_entities)
    entities["ro-crate-metadata.json"]["@type"] = "Thing"
    with pytest.raises(ValueError, match="must be of type"):
        find_root_entity_id(entities)
    # root type is not Dataset
    entities = deepcopy(orig_entities)
    entities["./"]["@type"] = "Thing"
    with pytest.raises(ValueError, match="must have"):
        find_root_entity_id(entities)


@pytest.mark.filterwarnings("ignore")
def test_find_root_multiple_entries():
    orig_entities = {
        "http://example.org/ro-crate-metadata.json": {
            "@id": "http://example.org/ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "http://example.org/"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
        },
        "http://example.org/": {
            "@id": "http://example.org/",
            "@type": "Dataset",
            "hasPart": [
                {"@id": "http://example.com/"},
                {"@id": "http://example.com/ro-crate-metadata.json"}
            ]
        },
        "http://example.com/ro-crate-metadata.json": {
            "@id": "http://example.com/ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "http://example.com/"},
            "conformsTo": {"@id": "https://w3id.com/ro/crate/1.1"},
        },
        "http://example.com/": {
            "@id": "http://example.com/",
            "@type": "Dataset",
        },
    }

    def check_finds_org(entities):
        m_id, r_id = find_root_entity_id(entities)
        assert m_id == "http://example.org/ro-crate-metadata.json"
        assert r_id == "http://example.org/"

    def check_picks_one(entities):
        m_id, r_id = find_root_entity_id(entities)
        assert m_id in [f"http://example.{_}/ro-crate-metadata.json" for _ in ("org", "com")]
        assert r_id in [f"http://example.{_}/" for _ in ("org", "com")]

    check_finds_org(orig_entities)
    # no root candidate contains the other one
    mod_entities = deepcopy(orig_entities)
    del mod_entities["http://example.org/"]["hasPart"]
    check_picks_one(mod_entities)
    # each root candidate contains the other one
    mod_entities = deepcopy(orig_entities)
    mod_entities["http://example.com/"]["hasPart"] = [
        {"@id": "http://example.org/"},
        {"@id": "http://example.org/ro-crate-metadata.json"}
    ]
    check_picks_one(mod_entities)
    # "about" does not reference the root entity
    mod_entities = deepcopy(orig_entities)
    for about in "http://google.com", {"@id": "http://google.com"}:
        mod_entities["http://example.com/ro-crate-metadata.json"]["about"] = about
        check_finds_org(mod_entities)
    # metadata type is not CreativeWork
    mod_entities = deepcopy(orig_entities)
    mod_entities["http://example.com/ro-crate-metadata.json"]["@type"] = "Thing"
    check_finds_org(mod_entities)
    # root type is not Dataset
    mod_entities = deepcopy(orig_entities)
    mod_entities["http://example.com/"]["@type"] = "Thing"
    check_finds_org(mod_entities)


def test_find_root_multiple_types():
    entities = {_["@id"]: _ for _ in [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
        },
        {
            "@id": "./",
            "@type": ["Dataset", "RepositoryCollection"],
        },
    ]}
    m_id, r_id = find_root_entity_id(entities)
    assert m_id == "ro-crate-metadata.json"
    assert r_id == "./"
    # "Dataset" not included
    del entities["./"]["@type"][0]
    with pytest.raises(ValueError):
        find_root_entity_id(entities)
    # Check we're not trying to be too clever
    entities["./"]["@type"] = "NotADataset"
    with pytest.raises(ValueError):
        find_root_entity_id(entities)
