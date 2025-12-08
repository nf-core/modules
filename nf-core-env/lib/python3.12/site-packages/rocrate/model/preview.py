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

import os
from pathlib import Path

from jinja2 import Template
from .file import File


class Preview(File):
    """
    RO-Crate preview file

    This object holds a preview of an RO Crate in HTML format_
    """
    BASENAME = "ro-crate-preview.html"

    def __init__(self, crate, source=None, properties=None):
        super().__init__(crate, source, self.BASENAME, properties=properties)

    def _empty(self):
        # default properties of the metadata entry
        val = {
            "@id": self.BASENAME,
            "@type": "CreativeWork",
            "about": {"@id": "./"}
        }
        return val

    def generate_html(self):
        base_path = os.path.abspath(os.path.dirname(__file__))
        template = open(
            os.path.join(base_path, '..', 'templates', 'preview_template.html.j2'),
            'r', encoding='utf-8'
        )
        src = Template(template.read())

        def template_function(func):
            src.globals[func.__name__] = func
            return func

        @template_function
        def stringify(a):
            if type(a) is list:
                return ', '.join(a)
            elif type(a) is str:
                return a
            else:
                if a._jsonld and a._jsonld['name']:
                    return a._jsonld['name']
                else:
                    return str(a)

        @template_function
        def is_object_list(a):
            if type(a) is list:
                for obj in a:
                    if obj is not str:
                        return True
            return False

        template.close()
        context_entities = []
        data_entities = []
        for entity in self.crate.contextual_entities:
            context_entities.append(entity._jsonld)
        for entity in self.crate.data_entities:
            data_entities.append(entity._jsonld)
        out_html = src.render(crate=self.crate, context=context_entities, data=data_entities)
        return out_html

    def stream(self, chunk_size=8192):
        if self.source:
            yield from super().stream()
        else:
            yield self.id, str.encode(self.generate_html(), encoding='utf-8')

    def _has_writeable_stream(self):
        return True

    def write(self, dest_base):
        write_path = Path(dest_base) / self.id
        super()._write_from_stream(write_path)
