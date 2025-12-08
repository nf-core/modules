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
import warnings

from .model import Metadata, LegacyMetadata


def read_metadata(metadata_path):
    """\
    Read an RO-Crate metadata file.

    Return a tuple of two elements: the context; a dictionary that maps entity
    ids to the entities themselves.
    """
    if isinstance(metadata_path, dict):
        metadata = metadata_path
    else:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
    try:
        context = metadata['@context']
        graph = metadata['@graph']
    except KeyError:
        raise ValueError(f"{metadata_path} must have a @context and a @graph")
    return context, {_["@id"]: _ for _ in graph}


def _check_descriptor(descriptor, entities):
    if descriptor["@type"] != "CreativeWork":
        raise ValueError('metadata descriptor must be of type "CreativeWork"')
    try:
        root = entities[descriptor["about"]["@id"]]
    except (KeyError, TypeError):
        raise ValueError("metadata descriptor does not reference the root entity")
    if ("Dataset" not in root["@type"] if isinstance(root["@type"], list) else root["@type"] != "Dataset"):
        raise ValueError('root entity must have "Dataset" among its types')
    return descriptor["@id"], root["@id"]


def find_root_entity_id(entities):
    """\
    Find metadata file descriptor and root data entity.

    Expects as input a dictionary that maps JSON entity IDs to the entities
    themselves (like the second element returned by read_metadata).

    Return a tuple of the corresponding identifiers (descriptor, root).
    If the entities are not found, raise KeyError. If they are found,
    but they don't satisfy the required constraints, raise ValueError.

    In the general case, the metadata file descriptor id can be an
    absolute URI whose last path segment is "ro-crate-metadata.json[ld]".
    Since there can be more than one such id in the crate, we need to
    choose among the corresponding (descriptor, root) entity pairs. First, we
    exclude those that don't satisfy other constraints, such as the
    descriptor entity being of type CreativeWork, etc.; if this doesn't
    leave us with a single pair, we try to pick one with a
    heuristic. Suppose we are left with the (m1, r1) and (m2, r2) pairs:
    if r1 is the actual root of this crate, then m2 and r2 are regular
    files in it, and as such they must appear in r1's hasPart; r2,
    however, is not required to have a hasPart property listing other
    files. Thus, we look for a pair whose root entity "contains" all
    descriptor entities from other pairs. If there is no such pair, or there
    is more than one, we just return an arbitrary pair.

    """
    descriptor = entities.get(Metadata.BASENAME, entities.get(LegacyMetadata.BASENAME))
    if descriptor:
        return _check_descriptor(descriptor, entities)
    candidates = []
    for id_, e in entities.items():
        basename = id_.rsplit("/", 1)[-1]
        if basename == Metadata.BASENAME or basename == LegacyMetadata.BASENAME:
            try:
                candidates.append(_check_descriptor(e, entities))
            except ValueError:
                pass
    if not candidates:
        raise KeyError("Metadata file descriptor not found")
    elif len(candidates) == 1:
        return candidates[0]
    else:
        warnings.warn("Multiple metadata file descriptors, will pick one with a heuristic")
        descriptor_ids = set(_[0] for _ in candidates)
        for m_id, r_id in candidates:
            try:
                root = entities[r_id]
                part_ids = set(_["@id"] for _ in root["hasPart"])
            except KeyError:
                continue
            if part_ids >= descriptor_ids - {m_id}:
                # if True for more than one candidate, this pick is arbitrary
                return m_id, r_id
        return candidates[0]  # fall back to arbitrary pick
