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

from io import RawIOBase


class MemoryBuffer(RawIOBase):
    """
    A buffer class that supports reading and writing binary data.
    The buffer automatically resets upon reading to make sure all data is read only once.
    """

    def __init__(self):
        self._buffer = b''

    def write(self, data):
        if self.closed:
            raise ValueError('write to closed file')
        self._buffer += data
        return len(data)

    def read(self, size=-1):
        if self.closed:
            raise ValueError('read from closed file')
        if size < 0:
            data = self._buffer
            self._buffer = b''
        else:
            data = self._buffer[:size]
            self._buffer = self._buffer[size:]
        return data

    def __len__(self):
        return len(self._buffer)
