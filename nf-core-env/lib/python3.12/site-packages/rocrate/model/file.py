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

from pathlib import Path
import requests
import shutil
import urllib.request
import warnings
from io import BytesIO, StringIO
from urllib.parse import unquote

from .file_or_dir import FileOrDir
from ..utils import is_url, iso_now, Mode


class File(FileOrDir):

    def _empty(self):
        val = {
            "@id": self.id,
            "@type": 'File'
        }
        return val

    def _has_writeable_stream(self):
        if isinstance(self.source, (BytesIO, StringIO)):
            return True
        elif is_url(str(self.source)):
            return self.fetch_remote
        else:
            return self.source is not None

    def _write_from_stream(self, out_file_path):
        if not self._has_writeable_stream():
            # is this does not correspond to a writeable stream (i.e. it is a url but fetch_remote is False),
            # we still want to consume the stream to consume file headers, run the size calculation, etc.
            all(self.stream())
            return

        out_file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_file_path, 'wb') as out_file:
            for _, chunk in self.stream():
                out_file.write(chunk)

    def _copy_file(self, path, out_file_path):
        out_file_path.parent.mkdir(parents=True, exist_ok=True)
        if not out_file_path.exists() or not out_file_path.samefile(path):
            shutil.copy(path, out_file_path)
        if self.record_size:
            self._jsonld['contentSize'] = str(out_file_path.stat().st_size)

    def write(self, base_path):
        out_file_path = Path(base_path) / unquote(self.id)
        if isinstance(self.source, (BytesIO, StringIO)) or is_url(str(self.source)):
            self._write_from_stream(out_file_path)
        elif self.source is None:
            # Allows to record a File entity whose @id does not exist, see #73
            warnings.warn(f"No source for {self.id}")
        else:
            if self.crate.mode == Mode.READ:
                in_file_path = unquote(str(self.source))
            else:
                in_file_path = self.source
            self._copy_file(in_file_path, out_file_path)

    def _stream_from_stream(self, stream):
        size = 0
        read = stream.read()
        if isinstance(self.source, StringIO):
            read = read.encode('utf-8')
        while len(read) > 0:
            yield self.id, read
            size += len(read)
            read = stream.read()
            if isinstance(self.source, StringIO):
                read = read.encode('utf-8')

        if self.record_size:
            self._jsonld['contentSize'] = str(size)

    def _stream_from_url(self, url, chunk_size=8192):
        if self.fetch_remote or self.validate_url:
            if self.validate_url:
                if url.startswith("http"):
                    with requests.head(url) as response:
                        self._jsonld.update({
                            'contentSize': response.headers.get('Content-Length'),
                            'encodingFormat': response.headers.get('Content-Type')
                        })
                    if not self.fetch_remote:
                        date_published = response.headers.get("Last-Modified", iso_now())
                        self._jsonld['sdDatePublished'] = date_published
            if self.fetch_remote:
                size = 0
                self._jsonld['contentUrl'] = str(url)
                with urllib.request.urlopen(url) as response:
                    while chunk := response.read(chunk_size):
                        yield self.id, chunk
                        size += len(chunk)

                # yield once for an empty file
                if size == 0:
                    yield self.id, b""

                if self.record_size:
                    self._jsonld['contentSize'] = str(size)

    def _stream_from_file(self, path, chunk_size=8192):
        size = 0
        with open(path, 'rb') as f:
            while chunk := f.read(chunk_size):
                yield unquote(self.id), chunk
                size += len(chunk)

        # yield once for an empty file
        if size == 0:
            yield self.id, b""

        if self.record_size:
            self._jsonld['contentSize'] = str(size)

    def stream(self, chunk_size=8192):
        if isinstance(self.source, (BytesIO, StringIO)):
            yield from self._stream_from_stream(self.source)
        elif is_url(str(self.source)):
            yield from self._stream_from_url(self.source, chunk_size)
        elif self.source is None:
            # Allows to record a File entity whose @id does not exist, see #73
            warnings.warn(f"No source for {self.id}")
        else:
            if self.crate.mode == Mode.READ:
                path = unquote(str(self.source))
            else:
                path = self.source
            yield from self._stream_from_file(path, chunk_size)
