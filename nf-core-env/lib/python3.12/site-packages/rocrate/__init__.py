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

"""
Create/parse RO-Crate metadata.

This module intends to help create or parse
RO-Crate metadata, see rocrate_

.. _rocrate: https://w3id.org/ro/crate/
"""

__author__ = ", ".join((
    'Daniel Bauer',
    'Eli Chadwick',
    'Paul De Geest',
    'Bert Droesbeke',
    'Ignacio Eguinoa',
    'Alban Gaignard',
    'Matthias Hörtenhuber',
    'Sebastiaan Huber',
    'Bruno Kinoshita',
    'Simone Leo',
    'Luca Pireddu',
    'Laura Rodríguez-Navas',
    'Raül Sirvent',
    'Stian Soiland-Reyes',
    'Laurent Thomas'
))
__copyright__ = """\
Copyright 2019-2025 The University of Manchester, UK
Copyright 2020-2025 Vlaams Instituut voor Biotechnologie (VIB), BE
Copyright 2020-2025 Barcelona Supercomputing Center (BSC), ES
Copyright 2020-2025 Center for Advanced Studies, Research and Development in Sardinia (CRS4), IT
Copyright 2022-2025 École Polytechnique Fédérale de Lausanne, CH
Copyright 2024-2025 Data Centre, SciLifeLab, SE
Copyright 2024-2025 National Institute of Informatics (NII), JP
Copyright 2025 Senckenberg Society for Nature Research (SGN), DE
Copyright 2025 European Molecular Biology Laboratory (EMBL), Heidelberg, DE
"""
__license__ = ("Apache License, version 2.0 "
               "<https://www.apache.org/licenses/LICENSE-2.0>")

# for arcp scheme registration with urllib.parse
import arcp  # noqa

# Convenience export of public functions/types
from .model.metadata import Metadata  # noqa
from ._version import __version__  # noqa
