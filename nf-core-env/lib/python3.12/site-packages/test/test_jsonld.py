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

"""Tests for the JSON-LD methods of the ROCrate object."""

from uuid import uuid4

import pytest

from rocrate.rocrate import ROCrate


# --- add


def test_add_jsonld_raises_json_is_none():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_jsonld(None)


def test_add_jsonld_raises_json_is_empty():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_jsonld({})


def test_add_jsonld_raises_json_missing_id():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_jsonld({
            '@type': 'CreativeWork'
        })


def test_add_jsonld_raises_json_missing_type(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    with pytest.raises(ValueError, match='.*@type.*'):
        crate.add_jsonld({
            '@id': './',
            'license': 'NA'
        })


def test_add_jsonld_raises_json_duplicate_id(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    with pytest.raises(ValueError, match='.*already exists.*'):
        crate.add_jsonld({
            '@id': './',
            '@type': 'Dataset',
            'license': 'NA'
        })


def test_add_jsonld(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    new_entity_id = f'#{uuid4()}'

    crate.add_jsonld({
        '@id': new_entity_id,
        '@type': 'CreativeWork',
        'name': 'A test entity'
    })
    new_entity = crate.get(new_entity_id)
    assert new_entity.type == 'CreativeWork'
    assert new_entity["name"] == 'A test entity'


# --- update


def test_update_jsonld_raises_json_is_none():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.update_jsonld(None)


def test_update_jsonld_raises_json_is_empty():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.update_jsonld({})


def test_update_jsonld_raises_json_missing_id():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.update_jsonld({
            '@type': 'CreativeWork'
        })


def test_update_jsonld_raises_id_not_found(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    missing_entity_id = f'#{uuid4()}'

    with pytest.raises(ValueError, match='.*does not exist.*'):
        crate.update_jsonld({
            '@id': missing_entity_id,
            '@type': 'CreativeWork',
            'name': 'This entity does not exist in the RO-Crate'
        })


def test_update_jsonld(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    new_entity_id = f'#{uuid4()}'

    crate.add_jsonld({
        '@id': new_entity_id,
        '@type': 'CreativeWork',
        'name': 'A test entity'
    })

    entity_added = crate.get(new_entity_id)
    assert entity_added.type == 'CreativeWork'
    assert entity_added['name'] == 'A test entity'

    entity_added['name'] = 'No potatoes today'
    # N.B.: Properties that start with @ are ignored when updating.
    update_dict = entity_added._jsonld.copy()
    update_dict['@type'] = 'Dataset'
    crate.update_jsonld(update_dict)

    updated_entity = crate.get(new_entity_id)
    assert entity_added.id == updated_entity.id
    assert updated_entity.type == 'CreativeWork'
    assert updated_entity['name'] == 'No potatoes today'

    assert '@type' in update_dict


# --- add or update


def test_add_or_update_jsonld_raises_json_is_none():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_or_update_jsonld(None)


def test_add_or_update_jsonld_raises_json_is_empty():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_or_update_jsonld({})


def test_add_or_update_jsonld_raises_json_missing_id():
    crate = ROCrate()

    with pytest.raises(ValueError, match='.*non-empty JSON-LD.*'):
        crate.add_or_update_jsonld({
            '@type': 'CreativeWork'
        })


def test_add_or_update_add_jsonld(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    new_entity_id = f'#{uuid4()}'

    crate.add_or_update_jsonld({
        '@id': new_entity_id,
        '@type': 'CreativeWork',
        'name': 'A test entity'
    })
    new_entity = crate.get(new_entity_id)
    assert new_entity.type == 'CreativeWork'
    assert new_entity["name"] == 'A test entity'


def test_add_or_update_update_jsonld(test_data_dir):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    new_entity_id = f'#{uuid4()}'

    crate.add_jsonld({
        '@id': new_entity_id,
        '@type': 'CreativeWork',
        'name': 'A test entity'
    })

    entity_added = crate.get(new_entity_id)
    assert entity_added.type == 'CreativeWork'
    assert entity_added['name'] == 'A test entity'

    entity_added['name'] = 'No potatoes today'
    crate.add_or_update_jsonld(entity_added._jsonld.copy())

    updated_entity = crate.get(new_entity_id)
    assert entity_added.id == updated_entity.id
    assert entity_added.type == updated_entity.type
    assert updated_entity['name'] == 'No potatoes today'
