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

import io
import pytest
import requests
import os
import uuid
import sys
import zipfile
from itertools import product
from urllib.error import URLError

from rocrate.model import Dataset, Person
from rocrate.rocrate import ROCrate


@pytest.mark.parametrize("gen_preview,to_zip", [(False, False), (False, True), (True, False), (True, True)])
def test_file_writing(test_data_dir, tmpdir, helpers, gen_preview, to_zip):
    crate = ROCrate(gen_preview=gen_preview)
    crate_name = 'Test crate'
    crate.name = crate_name
    creator_id = '#001'
    creator_name = 'Lee Ritenour'
    new_person = Person(crate, creator_id, {'name': creator_name})
    crate.add(new_person)
    crate.creator = new_person

    sample_file_id = 'sample_file.txt'
    sample_file2_id = 'a/b/sample_file2.csv'
    test_dir_id = 'test_add_dir/'
    data_entity_ids = [sample_file_id, sample_file2_id, test_dir_id]
    file_subdir_id = 'sample_file_subdir.txt'

    sample_file = test_data_dir / sample_file_id
    file_returned = crate.add_file(sample_file, record_size=True)
    assert file_returned.id == sample_file_id
    file_returned_subdir = crate.add_file(sample_file, sample_file2_id)
    assert file_returned_subdir.id == sample_file2_id
    test_dir_path = test_data_dir / test_dir_id
    test_dir_entity = crate.add_directory(test_dir_path, test_dir_id)
    assert isinstance(test_dir_entity, Dataset)

    out_path = tmpdir / 'ro_crate_out'
    if to_zip:
        zip_path = tmpdir / 'ro_crate_out.crate.zip'
        crate.write_zip(zip_path)
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(out_path)
    else:
        out_path.mkdir()
        crate.write(out_path)

    metadata_path = out_path / helpers.METADATA_FILE_NAME
    assert metadata_path.exists()
    preview_path = out_path / helpers.PREVIEW_FILE_NAME
    if gen_preview:
        assert preview_path.exists()
    else:
        assert not preview_path.exists()
    file1 = out_path / sample_file_id
    file2 = out_path / sample_file2_id
    file_subdir = out_path / test_dir_id / file_subdir_id
    assert file1.exists()
    with open(sample_file) as f1, open(file1) as f2:
        sample_file_content = f1.read()
        assert sample_file_content == f2.read()
    assert file2.exists()
    with open(file2) as f:
        assert sample_file_content == f.read()
    assert file_subdir.exists()
    with open(test_dir_path / file_subdir_id) as f1, open(file_subdir) as f2:
        assert f1.read() == f2.read()

    json_entities = helpers.read_json_entities(out_path)
    helpers.check_crate(json_entities, data_entity_ids=data_entity_ids)
    root = json_entities["./"]
    assert root["name"] == crate_name
    assert "datePublished" in root
    assert root["creator"] == {"@id": creator_id}
    assert creator_id in json_entities
    assert json_entities[creator_id]["name"] == creator_name
    if gen_preview:
        assert helpers.PREVIEW_FILE_NAME in json_entities
    file_entity = json_entities[sample_file_id]
    assert file_entity["contentSize"] == str(file1.stat().st_size)


@pytest.mark.parametrize("stream_cls", [io.BytesIO, io.StringIO])
def test_in_mem_stream(stream_cls, tmpdir, helpers):
    crate = ROCrate()

    test_file_id = 'a/b/test_file.txt'
    file_content = b'\x00\x01\x02' if stream_cls is io.BytesIO else 'foo'
    file_returned = crate.add_file(stream_cls(file_content), test_file_id, record_size=True)
    assert file_returned.id == test_file_id

    out_path = tmpdir / 'ro_crate_out'
    out_path.mkdir()
    crate.write(out_path)

    metadata_path = out_path / helpers.METADATA_FILE_NAME
    assert metadata_path.exists()
    file1 = out_path / test_file_id
    assert file1.exists()
    mode = 'r' + ('b' if stream_cls is io.BytesIO else 't')
    with open(file1, mode) as f:
        assert f.read() == file_content
    json_entities = helpers.read_json_entities(out_path)
    file_entity = json_entities[test_file_id]
    assert file_entity['contentSize'] == '3'


@pytest.mark.parametrize(
    "fetch_remote,validate_url,to_zip",
    list(product((False, True), repeat=3))
)
def test_remote_uri(tmpdir, helpers, fetch_remote, validate_url, to_zip):
    crate = ROCrate()
    url = ('https://raw.githubusercontent.com/ResearchObject/ro-crate-py/'
           'master/test/test-data/sample_file.txt')
    relpath = "a/b/sample_file.txt"
    kw = {
        "fetch_remote": fetch_remote,
        "validate_url": validate_url,
        "record_size": True,
    }
    if fetch_remote:
        file_ = crate.add_file(url, relpath, **kw)
        assert file_.id == relpath
    else:
        file_ = crate.add_file(url, **kw)
        assert file_.id == url

    out_path = tmpdir / 'ro_crate_out'
    if to_zip:
        zip_path = f"{out_path}.zip"
        crate.write_zip(zip_path)
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(out_path)
    else:
        crate.write(out_path)

    out_crate = ROCrate(out_path)
    if fetch_remote:
        out_file = out_crate.dereference(file_.id)
        assert (out_path / relpath).is_file()
        assert out_file["contentUrl"] == url
        assert out_file["contentSize"] == str((out_path / file_.id).stat().st_size)
    else:
        out_file = out_crate.dereference(url)
        assert not (out_path / relpath).exists()
    if validate_url:
        props = out_file.properties()
        assert {"contentSize", "encodingFormat"}.issubset(props)
        if not fetch_remote:
            assert "sdDatePublished" in props


def test_file_uri(tmpdir):
    f_name = uuid.uuid4().hex
    f_path = (tmpdir / f_name).resolve()
    f_uri = f"file:///{f_path}"  # extra slash needed on some windows systems
    with open(f_path, "wt") as f:
        f.write("FOO\n")
    crate = ROCrate()
    f_entity = crate.add_file(f_uri, fetch_remote=True)
    assert f_entity.id == f_name

    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)

    out_crate = ROCrate(out_path)
    assert out_crate.dereference(f_name) is not None
    assert (out_path / f_name).is_file()


@pytest.mark.skipif(os.name != "posix", reason="':' not allowed in dir name")
def test_looks_like_file_uri(tmpdir, monkeypatch):
    f_name = uuid.uuid4().hex
    f_parent = (tmpdir / "file:")
    f_parent.mkdir()
    f_path = f_parent / f_name
    with open(f_path, "wt") as f:
        f.write("FOO\n")
    monkeypatch.chdir(tmpdir)
    crate = ROCrate()
    # Missing if interpreted as URI, present if intepreted as path
    uri = f"file:/{f_name}"
    entity = crate.add_file(uri, fetch_remote=False)
    assert entity.id == uri

    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)

    out_crate = ROCrate(out_path)
    assert out_crate.dereference(uri) is not None
    assert not (out_path / f_name).is_file()
    assert not (out_path / uri).is_file()

    # Check that the file can be added to the crate using its absolute path
    entity = crate.add_file(f_path.resolve())
    assert entity.id == f_name

    out_path = tmpdir / 'ro_crate_out_updated'
    crate.write(out_path)

    out_crate = ROCrate(out_path)
    assert out_crate.dereference(f_name) is not None
    assert (out_path / f_name).is_file()


@pytest.mark.skipif(sys.platform == "darwin", reason="FTP opening broken on macOS CI instances")
@pytest.mark.slow
@pytest.mark.parametrize("fetch_remote", [False, True])
def test_ftp_uri(tmpdir, fetch_remote):
    crate = ROCrate()
    url = 'ftp://ftp-trace.ncbi.nih.gov/pub/gdp/README'
    relpath = "a/b/README"
    if fetch_remote:
        file_ = crate.add_file(url, relpath, fetch_remote=True)
        assert file_.id == relpath
    else:
        file_ = crate.add_file(url, fetch_remote=False)
        assert file_.id == url

    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)
    if fetch_remote:
        assert (out_path / relpath).is_file()
    else:
        assert not (out_path / relpath).exists()


@pytest.mark.parametrize(
    "fetch_remote,validate_url",
    list(product((False, True), repeat=2))
)
def test_remote_dir(tmpdir, helpers, fetch_remote, validate_url):
    crate = ROCrate()
    url = "https://ftp.mozilla.org/pub/misc/errorpages/"
    relpath = "pub/misc/errorpages/"
    properties = {
        "hasPart": [
            {"@id": "404.html"},
            {"@id": "500.html"},
        ],
    }
    kw = {
        "fetch_remote": fetch_remote,
        "validate_url": validate_url,
        "properties": properties,
    }
    if fetch_remote:
        dataset = crate.add_dataset(url, relpath, **kw)
        assert dataset.id == relpath
    else:
        dataset = crate.add_dataset(url, **kw)
        assert dataset.id == url

    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)

    out_crate = ROCrate(out_path)
    if fetch_remote:
        out_dataset = out_crate.dereference(relpath)
        assert (out_path / relpath).is_dir()
        for entry in properties["hasPart"]:
            assert (out_path / relpath / entry["@id"]).is_file()
    else:
        out_dataset = out_crate.dereference(url)
        assert not (out_path / relpath).exists()
    if validate_url and not fetch_remote:
        assert "sdDatePublished" in out_dataset.properties()


def test_remote_uri_exceptions(tmpdir):
    crate = ROCrate()
    url = ('https://raw.githubusercontent.com/ResearchObject/ro-crate-py/'
           f'master/test/test-data/{uuid.uuid4().hex}.foo')
    crate.add_file(source=url, fetch_remote=True)
    out_path = tmpdir / 'ro_crate_out_1'
    out_path.mkdir()
    with pytest.raises(URLError):
        crate.write(out_path)

    crate = ROCrate()
    url = ('https://raw.githubusercontent.com/ResearchObject/ro-crate-py/'
           'master/test/test-data/sample_file.txt')
    crate.add_file(source=url, dest_path="a/sample_file.txt", fetch_remote=True)
    out_path = tmpdir / 'ro_crate_out_2'
    out_path.mkdir()
    (out_path / "a").mkdir(mode=0o444)
    try:
        crate.write(out_path)
    except PermissionError:
        pass
    # no error on Windows, or on Linux as root, so we don't use pytest.raises


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize("what", ["file", "dataset"])
def test_missing_source(test_data_dir, tmpdir, what):
    path = test_data_dir / uuid.uuid4().hex

    crate = ROCrate()
    entity = getattr(crate, f"add_{what}")(path)
    assert entity is crate.dereference(path.name)
    crate_dir = tmpdir / 'ro_crate_out_1'
    with pytest.raises(OSError):
        crate.write(crate_dir)

    crate = ROCrate()
    entity = getattr(crate, f"add_{what}")(path, path.name)
    assert entity is crate.dereference(path.name)
    crate_dir = tmpdir / 'ro_crate_out_2'
    with pytest.raises(OSError):
        crate.write(crate_dir)

    crate = ROCrate()
    entity = getattr(crate, f"add_{what}")(None, path.name)
    assert entity is crate.dereference(path.name)
    crate_dir = tmpdir / 'ro_crate_out_3'
    crate.write(crate_dir)
    out_path = crate_dir / path.name
    if what == "file":
        assert not out_path.exists()
    else:
        assert out_path.is_dir()
        assert not any(out_path.iterdir())


@pytest.mark.parametrize("fetch_remote,validate_url", [(False, False), (False, True), (True, False), (True, True)])
def test_stringio_no_dest(test_data_dir, fetch_remote, validate_url):
    crate = ROCrate()
    with pytest.raises(ValueError):
        crate.add_file(io.StringIO("foo"))


@pytest.mark.parametrize("fetch_remote,validate_url", [(False, False), (False, True), (True, False), (True, True)])
def test_no_source_no_dest(test_data_dir, fetch_remote, validate_url):
    crate = ROCrate()
    with pytest.raises(ValueError):
        crate.add_file()


def test_dataset(test_data_dir, tmpdir):
    crate = ROCrate()
    path_a_b = test_data_dir / "a" / "b"
    path_c = test_data_dir / "c"
    for p in path_a_b, path_c:
        p.mkdir(parents=True)
    d1 = crate.add_dataset(path_a_b)
    assert crate.dereference("b") is d1
    d2 = crate.add_dataset(path_a_b, "a/b")
    assert crate.dereference("a/b") is d2
    d_from_str = crate.add_dataset(str(path_c))
    assert crate.dereference("c") is d_from_str

    out_path = tmpdir / 'ro_crate_out'
    out_path.mkdir()
    crate.write(out_path)

    assert (out_path / "b").is_dir()
    assert (out_path / "a" / "b").is_dir()
    assert (out_path / "c").is_dir()


def test_no_parts(tmpdir, helpers):
    crate = ROCrate()

    out_path = tmpdir / 'ro_crate_out'
    out_path.mkdir()
    crate.write(out_path)

    json_entities = helpers.read_json_entities(out_path)
    helpers.check_crate(json_entities)
    assert "hasPart" not in json_entities["./"]


def test_write_zip_copy_unlisted(test_data_dir, tmpdir):
    crate_dir = test_data_dir / 'ro-crate-galaxy-sortchangecase'
    crate = ROCrate(crate_dir)

    zip_name = 'ro_crate_out.crate.zip'
    zip_path = tmpdir / zip_name
    crate.write_zip(zip_path)
    out_path = tmpdir / 'ro_crate_out'
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(out_path)

    assert (out_path / "test" / "test1" / "input.bed").is_file()
    assert (out_path / "test" / "test1" / "output_exp.bed").is_file()


def test_no_zip_in_zip(test_data_dir, tmpdir):
    crate_dir = test_data_dir / 'ro-crate-galaxy-sortchangecase'
    crate = ROCrate(crate_dir)

    zip_name = 'ro_crate_out.crate.zip'
    zip_path = crate_dir / zip_name  # within the crate dir
    crate.write_zip(zip_path)
    out_path = tmpdir / 'ro_crate_out'
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(out_path)

    assert not (out_path / zip_name).exists()


def test_add_tree(test_data_dir, tmpdir):
    source = test_data_dir / "read_extra"
    (source / "ro-crate-metadata.json").unlink()
    crate = ROCrate()
    crate.add_tree(source)
    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)

    assert (out_path / "read_extra").is_dir()
    assert (out_path / "read_extra" / "listed").is_dir()
    assert (out_path / "read_extra" / "listed" / "listed.txt").is_file()
    assert (out_path / "read_extra" / "listed" / "not_listed.txt").is_file()
    assert (out_path / "read_extra" / "not_listed").is_dir()
    assert (out_path / "read_extra" / "not_listed" / "not_listed.txt").is_file()
    assert (out_path / "read_extra" / "listed.txt").is_file()
    assert (out_path / "read_extra" / "not_listed.txt").is_file()

    crate = ROCrate(out_path)
    top = crate.get("read_extra")
    assert top.type == "Dataset"
    listed = crate.get("read_extra/listed")
    assert listed.type == "Dataset"
    not_listed = crate.get("read_extra/not_listed")
    assert not_listed.type == "Dataset"
    listed_txt = crate.get("read_extra/listed.txt")
    assert listed_txt.type == "File"
    not_listed_txt = crate.get("read_extra/not_listed.txt")
    assert not_listed_txt.type == "File"
    assert set(top["hasPart"]) == {listed, not_listed, listed_txt, not_listed_txt}
    listed_listed_txt = crate.get("read_extra/listed/listed.txt")
    listed_not_listed_txt = crate.get("read_extra/listed/not_listed.txt")
    assert set(listed["hasPart"]) == {listed_listed_txt, listed_not_listed_txt}
    not_listed_not_listed_txt = crate.get("read_extra/not_listed/not_listed.txt")
    assert set(not_listed["hasPart"]) == {not_listed_not_listed_txt}

    with pytest.raises(ValueError):
        crate.add_tree(None, dest_path="foobar")


def test_http_header(tmpdir):
    crate = ROCrate()
    url = "https://zenodo.org/records/10782431/files/lysozyme_datasets.zip"
    file_ = crate.add_file(url, validate_url=True)
    assert file_.id == url
    out_path = tmpdir / 'ro_crate_out'
    crate.write(out_path)
    out_crate = ROCrate(out_path)
    out_file = out_crate.dereference(url)
    props = out_file.properties()
    assert props.get("encodingFormat") == "application/octet-stream"
    assert "sdDatePublished" in props
    with requests.head(url) as response:
        assert props["sdDatePublished"] == response.headers.get("last-modified")


def test_stream(test_data_dir, tmpdir):
    source = test_data_dir / "read_crate"
    crate = ROCrate(source)

    out_path = tmpdir / 'ro_crate_out.zip'
    with open(out_path, "wb") as out:
        for chunk in crate.stream_zip():
            out.write(chunk)

    with zipfile.ZipFile(out_path, "r") as zf:
        assert not zf.testzip()
        for info in zf.infolist():
            assert info.file_size > 0

    extract_path = tmpdir / 'ro_crate_out'
    with zipfile.ZipFile(out_path, "r") as zf:
        zf.extractall(extract_path)
    assert (extract_path / "ro-crate-metadata.jsonld").is_file()
    assert (extract_path / "examples" / "README.txt").is_file()
    assert (extract_path / "test" / "test-metadata.json").is_file()


def test_percent_escape(test_data_dir, tmpdir, helpers):
    crate = ROCrate()
    f_path = test_data_dir / "read_crate" / "with space.txt"
    f1 = crate.add_file(f_path)
    assert f1.id == "with%20space.txt"
    f2 = crate.add_file(f_path, dest_path="subdir/with space.txt")
    assert f2.id == "subdir/with%20space.txt"
    f3 = crate.add_file(test_data_dir / "read_crate" / "without%20space.txt")
    assert f3.id == "without%2520space.txt"
    d_path = test_data_dir / "read_crate" / "a b"
    d1 = crate.add_dataset(d_path)
    assert d1.id == "a%20b/"
    d2 = crate.add_dataset(d_path, dest_path="subdir/a b")
    assert d2.id == "subdir/a%20b/"
    d3 = crate.add_dataset(test_data_dir / "read_crate" / "j%20k")
    assert d3.id == "j%2520k/"
    out_path = tmpdir / "ro_crate_out"
    crate.write(out_path)
    json_entities = helpers.read_json_entities(out_path)
    assert "with%20space.txt" in json_entities
    assert "subdir/with%20space.txt" in json_entities
    assert "without%2520space.txt" in json_entities
    assert "a%20b/" in json_entities
    assert "subdir/a%20b/" in json_entities
    assert "j%2520k/" in json_entities
    assert (out_path / "with space.txt").is_file()
    assert (out_path / "subdir" / "with space.txt").is_file()
    assert (out_path / "without%20space.txt").is_file()
    assert (out_path / "a b" / "c d.txt").is_file()
    assert (out_path / "subdir" / "a b" / "c d.txt").is_file()
    assert (out_path / "j%20k" / "l%20m.txt").is_file()
    out_zip_path = tmpdir / "ro_crate_out.zip"
    crate.write_zip(out_zip_path)
    unpack_path = tmpdir / "unpack"
    with zipfile.ZipFile(out_zip_path, "r") as zf:
        zf.extractall(unpack_path)
    json_entities = helpers.read_json_entities(unpack_path)
    assert "with%20space.txt" in json_entities
    assert "subdir/with%20space.txt" in json_entities
    assert "without%2520space.txt" in json_entities
    assert "a%20b/" in json_entities
    assert "subdir/a%20b/" in json_entities
    assert "j%2520k/" in json_entities
    assert (unpack_path / "with space.txt").is_file()
    assert (unpack_path / "subdir" / "with space.txt").is_file()
    assert (unpack_path / "without%20space.txt").is_file()
    assert (unpack_path / "a b" / "c d.txt").is_file()
    assert (unpack_path / "subdir" / "a b" / "c d.txt").is_file()
    assert (unpack_path / "j%20k" / "l%20m.txt").is_file()


def test_stream_empty_file(test_data_dir, tmpdir):
    """
    Test that empty files are written correctly to the zip file.
    """
    crate = ROCrate()
    crate_dir = test_data_dir / "empty_file_crate"
    crate.add_file(crate_dir / "empty.txt")
    crate.add_directory(crate_dir / "folder")

    # write the crate to a zip file
    out_path = tmpdir / 'ro_crate_out.zip'
    crate.write_zip(out_path)

    # Check that the zip file contains empty files
    assert out_path.is_file()
    files_in_zip = {}
    with zipfile.ZipFile(out_path, "r") as zf:
        for info in zf.infolist():
            files_in_zip[info.filename] = info.file_size

    assert files_in_zip["empty.txt"] == 0
    assert files_in_zip["folder/empty_not_listed.txt"] == 0


def test_write_zip_nested_dest(tmpdir, helpers):
    root = tmpdir / "root"
    root.mkdir()
    (root / "a b").mkdir()
    (root / "a b" / "c d.txt").write_text("C D\n")
    (root / "a b" / "j k").mkdir()
    (root / "a b" / "j k" / "l m.txt").write_text("L M\n")
    crate = ROCrate()
    d1 = crate.add_dataset(root / "a b", dest_path="subdir/a b")
    assert d1.id == "subdir/a%20b/"
    out_zip_path = tmpdir / "ro_crate_out.zip"
    crate.write_zip(out_zip_path)
    unpack_path = tmpdir / "unpack"
    with zipfile.ZipFile(out_zip_path, "r") as zf:
        zf.extractall(unpack_path)
    json_entities = helpers.read_json_entities(unpack_path)
    assert "subdir/a%20b/" in json_entities
    assert (unpack_path / "subdir" / "a b" / "c d.txt").is_file()
    assert (unpack_path / "subdir" / "a b" / "j k" / "l m.txt").is_file()
