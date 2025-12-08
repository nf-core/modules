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


from .entity import Entity


class DataEntity(Entity):

    def write(self, base_path):
        pass

    def stream(self, chunk_size=8192):
        """ Stream the data from the source. Each chunk of the content is yielded as a tuple
        containing the name of the destination file relative to the crate and the chunk of data.
        The destination file name is required because a DataEntity can be a file or a
        collection of files (Dataset) and the caller need to know to which file a chunk belongs.
        For collection of files, the caller can assume that files are streamed one after another,
        meaning once the destination name changes, a file can be closed and the next one can be
        openend.
        """
        yield from ()
