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

import os
from pathlib import Path

from rocrate.rocrate import ROCrate


def test_crate_update(test_data_dir, tmpdir, helpers):
    crate_dir = test_data_dir / 'read_crate'
    crate = ROCrate(crate_dir)

    content_map = {}
    for root, dirs, files in os.walk(crate_dir):
        root = Path(root)
        for name in files:
            if not name.startswith("ro-crate"):
                path = root / name
                with open(path, "rb") as f:
                    content_map[path] = f.read()

    # update an existing file
    upd_file_id = "test_file_galaxy.txt"
    upd_file = crate.dereference(upd_file_id)
    with open(upd_file.source, "rb") as f:
        content = f.read()
    upd_content = content + b"foobar\n"
    upd_source = tmpdir / "upd_source.txt"
    with open(upd_source, "wb") as f:
        f.write(upd_content)
    crate.delete(upd_file)
    crate.add_file(upd_source, upd_file_id)

    # add a new file
    new_file_id = "spam.txt"
    new_content = b"enlarge your crate\n"
    new_source = tmpdir / "new_source.txt"
    with open(new_source, "wb") as f:
        f.write(new_content)
    crate.add_file(new_source, new_file_id)

    crate.write(crate_dir)

    for root, dirs, files in os.walk(crate_dir):
        root = Path(root)
        for name in files:
            if not name.startswith("ro-crate"):
                path = root / name
                with open(path, "rb") as f:
                    content = f.read()
                if path == crate_dir / upd_file_id:
                    assert content == upd_content
                elif path == crate_dir / new_file_id:
                    assert content == new_content
                else:
                    assert content == content_map[path]
