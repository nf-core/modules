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

"""
Model of RO-Crate metadata.

This module intends to cover each of the data entities and contextual entities
in rocrate_ represented as different Python classes.

.. _rocrate: https://w3id.org/ro/crate/
"""

from .computationalworkflow import ComputationalWorkflow, WorkflowDescription, Workflow
from .computerlanguage import ComputerLanguage
from .contextentity import ContextEntity
from .creativework import CreativeWork
from .data_entity import DataEntity
from .dataset import Dataset
from .entity import Entity
from .file import File
from .file_or_dir import FileOrDir
from .metadata import Metadata, LegacyMetadata
from .person import Person
from .root_dataset import RootDataset
from .softwareapplication import SoftwareApplication
from .testdefinition import TestDefinition
from .testinstance import TestInstance
from .preview import Preview
from .testservice import TestService
from .testsuite import TestSuite

__all__ = [
    "ComputationalWorkflow",
    "ComputerLanguage",
    "ContextEntity",
    "CreativeWork",
    "DataEntity",
    "Dataset",
    "Entity",
    "File",
    "FileOrDir",
    "LegacyMetadata",
    "Metadata",
    "Person",
    "Preview",
    "RootDataset",
    "SoftwareApplication",
    "TestDefinition",
    "TestInstance",
    "TestService",
    "TestSuite",
    "Workflow",
    "WorkflowDescription",
]
