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

from .contextentity import ContextEntity


class TestSuite(ContextEntity):

    def _empty(self):
        return {
            "@id": self.id,
            "@type": 'TestSuite'
        }

    @property
    def _default_type(self):
        return "TestSuite"

    @property
    def name(self):
        return self.get("name")

    @name.setter
    def name(self, name):
        self["name"] = name

    @property
    def instance(self):
        return self.get("instance")

    @instance.setter
    def instance(self, instance):
        self["instance"] = instance

    @property
    def definition(self):
        return self.get("definition")

    @definition.setter
    def definition(self, definition):
        self["definition"] = definition
