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

import datetime
import json
import sys
import tempfile
import timeit
import uuid
from pathlib import Path

import pytest
from rocrate.rocrate import ROCrate
from rocrate.model import (
    DataEntity,
    File,
    Dataset,
    ComputationalWorkflow,
    Person,
    Preview,
    ContextEntity
)


RAW_REPO_URL = "https://raw.githubusercontent.com/ResearchObject/ro-crate-py"


def test_dereferencing(test_data_dir, helpers, monkeypatch):
    crate = ROCrate(gen_preview=True)

    # verify default entities
    root_dataset = crate.dereference('./')
    assert crate.root_dataset is root_dataset
    metadata_entity = crate.dereference(helpers.METADATA_FILE_NAME)
    assert crate.metadata is metadata_entity
    preview_entity = crate.dereference(helpers.PREVIEW_FILE_NAME)
    assert preview_entity is crate.preview

    # dereference added files
    sample_file = test_data_dir / 'sample_file.txt'
    file_returned = crate.add_file(sample_file)
    assert isinstance(file_returned, File)
    dereference_file = crate.dereference("sample_file.txt")
    assert file_returned is dereference_file
    readme_url = f'{RAW_REPO_URL}/master/README.md'
    readme_entity = crate.add_file(readme_url)
    assert crate.dereference(readme_url) is readme_entity

    # dereference nonexistent entity
    assert crate.dereference(f"#{uuid.uuid4()}") is None
    assert crate.dereference(f"#{uuid.uuid4()}", file_returned) is file_returned

    # alias
    assert crate.get(readme_url) is readme_entity

    # dereference local dir added with trailing slash
    monkeypatch.chdir(test_data_dir)
    local_dir = "test_add_dir/"
    local_dir_entity = crate.add_dataset(local_dir)
    assert crate.dereference(local_dir) is local_dir_entity


@pytest.mark.parametrize("name", [".foo", "foo.", ".foo/", "foo./"])
def test_dereferencing_equivalent_id(test_data_dir, name):
    test_ids = [name, f'./{name}', f'.//{name}', f'./a/../{name}']
    path = test_data_dir / name
    if name.endswith("/"):
        path.mkdir()
        test_ids.extend([f"{name}/", f"{name}//"])
    else:
        path.touch()
    for id_ in test_ids:
        crate = ROCrate()
        if name.endswith("/"):
            entity = crate.add_directory(path, name)
        else:
            entity = crate.add_file(path, name)
        deref_entity = crate.dereference(id_)
        assert deref_entity is entity


def test_data_entities(test_data_dir):
    crate = ROCrate()
    file_ = crate.add(File(crate, test_data_dir / 'sample_file.txt'))
    dataset = crate.add(Dataset(crate, test_data_dir / 'test_add_dir'))
    data_entity = crate.add(DataEntity(crate, '#mysterious'))
    assert set(crate.data_entities) == {file_, dataset, data_entity}
    part_ids = set(_["@id"] for _ in crate.root_dataset._jsonld["hasPart"])
    assert set(_.id for _ in (file_, dataset, data_entity)) <= part_ids


@pytest.mark.skipif(sys.platform == "darwin", reason="CI sometimes fails on macOS")
def test_data_entities_perf():
    """\
    Test that adding a data entity happens in constant time.

    See https://github.com/ResearchObject/ro-crate-py/pull/127 (this is a
    regression test). The time required for 500 iterations should be ~0.01s if
    entities are added in constant time, ~1s if they are added in linear time.
    """
    crate = ROCrate()
    assert timeit.Timer(
        "crate.add_file(uuid.uuid4().hex)", "import uuid", globals=locals()
    ).timeit(500) < 0.1


def test_remote_data_entities():
    crate = ROCrate()
    file_uri = "https://www.rfc-editor.org/rfc/rfc3986.txt"
    dataset_uri = "https://ftp.mozilla.org/pub/misc/errorpages/"
    file_ = crate.add(File(crate, file_uri))
    dataset = crate.add(Dataset(crate, dataset_uri))
    assert file_.id == file_uri
    assert dataset.id == dataset_uri


def test_bad_data_entities(test_data_dir):
    # no source and no dest_path
    crate = ROCrate()
    with pytest.raises(ValueError):
        crate.add(Dataset(crate))
    with pytest.raises(ValueError):
        crate.add(File(crate))
    # absolute dest_path
    with pytest.raises(ValueError):
        tmp = Path(tempfile.gettempdir())
        crate.add(Dataset(crate, test_data_dir, tmp / "foo"))
    with pytest.raises(ValueError):
        crate.add(File(crate, test_data_dir / "sample_file.txt", tmp / "x.txt"))


def test_contextual_entities():
    crate = ROCrate()
    new_person = crate.add(Person(crate, '#joe', {'name': 'Joe Pesci'}))
    person_dereference = crate.dereference('#joe')
    assert person_dereference is new_person
    assert person_dereference.type == 'Person'
    person_prop = person_dereference.properties()
    assert person_prop['@type'] == 'Person'
    assert person_prop['name'] == 'Joe Pesci'
    assert not new_person.datePublished
    id_ = "https://orcid.org/0000-0002-1825-0097"
    new_person = crate.add(Person(crate, id_, {'name': 'Josiah Carberry'}))
    person_dereference = crate.dereference(id_)
    assert person_dereference is new_person


def test_properties():
    crate = ROCrate()

    crate_name = "new crate"
    crate.name = crate_name
    assert crate.name == crate_name

    crate_description = "this is a new crate"
    crate.description = crate_description
    assert crate.description == crate_description

    assert crate.datePublished == crate.root_dataset.datePublished
    assert isinstance(crate.root_dataset.datePublished, datetime.datetime)
    assert isinstance(crate.root_dataset["datePublished"], str)
    crate_datePublished = datetime.datetime.now()
    crate.datePublished = crate_datePublished
    assert crate.datePublished == crate_datePublished

    new_person = crate.add(Person(crate, '#001', {'name': 'Lee Ritenour'}))
    crate.creator = new_person
    assert crate.creator is new_person
    assert isinstance(crate.creator, Person)
    assert crate.creator['name'] == 'Lee Ritenour'
    assert crate.creator.type == 'Person'

    new_person2 = crate.add(Person(crate, '#002', {'name': 'Lee Ritenour'}))
    crate.creator = [new_person, new_person2]
    assert isinstance(crate.creator, list)
    assert crate.creator[0] is new_person
    assert crate.creator[1] is new_person2


def test_uuid():
    crate = ROCrate()
    new_person = crate.add(Person(crate, properties={"name": "No Identifier"}))
    jsonld = new_person.as_jsonld()
    assert "Person" == jsonld["@type"]
    assert jsonld["@id"].startswith("#")
    # Check it made a valid UUIDv4
    u = uuid.UUID(jsonld["@id"][1:])
    assert 4 == u.version


def test_update(test_data_dir, tmpdir, helpers):
    crate = ROCrate()
    assert "hasPart" not in crate.root_dataset
    wf_path = test_data_dir / "test_galaxy_wf.ga"
    john = crate.add(Person(crate, '#john', {'name': 'John Doe'}))
    file_ = crate.add_file(wf_path)
    assert crate.root_dataset["hasPart"] == [file_]
    file_["author"] = john
    assert isinstance(file_, File)
    assert crate.mainEntity is None
    wf = crate.add_workflow(wf_path, main=True, lang="galaxy")
    assert isinstance(wf, ComputationalWorkflow)
    assert "author" not in wf
    assert crate.mainEntity is wf
    assert crate.dereference(john.id) is john
    assert crate.dereference(file_.id) is wf
    assert crate.root_dataset["hasPart"] == [wf]
    assert crate.root_dataset.properties()["hasPart"] == [{"@id": "test_galaxy_wf.ga"}]

    out_path = tmpdir / "ro_crate_out"
    out_path.mkdir()
    crate.write(out_path)
    json_entities = helpers.read_json_entities(out_path)
    helpers.check_wf_crate(json_entities, wf.id)


def test_delete(test_data_dir):
    crate = ROCrate()
    # fundamental entities
    with pytest.raises(ValueError):
        crate.delete(crate.root_dataset)
    with pytest.raises(ValueError):
        crate.delete(crate.metadata)
    # preview
    preview = crate.add(Preview(crate))
    assert preview in crate.default_entities
    crate.delete(preview)
    assert preview not in crate.default_entities
    assert crate.preview is None
    # data entities
    file1 = crate.add_file(test_data_dir / "sample_file.txt")
    file2 = crate.add_file(test_data_dir / "sample_cwl_wf.cwl")
    for f in file1, file2:
        assert f in crate.root_dataset["hasPart"]
        assert f in crate.data_entities
    crate.delete(file1)
    assert file1 not in crate.data_entities
    assert file1 not in crate.root_dataset["hasPart"]
    crate.delete(file2)
    assert file2 not in crate.data_entities
    assert "hasPart" not in crate.root_dataset
    crate.delete(file2)  # no-op
    # contextual entities
    john = crate.add(Person(crate, '#john', {'name': 'John Doe'}))
    assert john in crate.contextual_entities
    crate.delete(john)
    assert john not in crate.contextual_entities
    crate.delete(john)  # no-op


# FIXME: what to do with refs is still WIP
def test_delete_refs(test_data_dir, tmpdir, helpers):
    def_path = "test/test1/sort-and-change-case-test.yml"
    crate = ROCrate(test_data_dir / 'ro-crate-galaxy-sortchangecase')
    suite = crate.dereference("#test1")
    definition = crate.dereference(def_path)
    assert suite.definition is definition
    crate.delete(definition)
    assert suite.definition is not definition  # so far, so good
    assert suite.definition == str(def_path)  # should probably be set to None
    # check json output
    out_path = tmpdir / "ro_crate_out"
    crate.write(out_path)
    json_entities = helpers.read_json_entities(out_path)
    assert def_path not in json_entities  # good
    assert json_entities["#test1"]["definition"]["@id"] == def_path  # not good
    # the test suite should not be in the crate at all
    # even better, the write should fail in such an inconsistent state
    # or perhaps such an inconsistent state should not be allowed at all


def test_delete_by_id(test_data_dir):
    crate = ROCrate()
    path = test_data_dir / "sample_file.txt"
    f = crate.add_file(path)  # with this signature, the file's id will be its basename
    assert f in crate.data_entities
    assert f in crate.root_dataset["hasPart"]
    assert f.id == path.name
    crate.delete(path.name)
    assert f not in crate.data_entities
    assert "hasPart" not in crate.root_dataset


def test_self_delete(test_data_dir):
    crate = ROCrate()
    path = test_data_dir / "sample_file.txt"
    f = crate.add_file(path)  # with this signature, the file's id will be its basename
    assert f in crate.data_entities
    assert f in crate.root_dataset["hasPart"]
    f.delete()
    assert f not in crate.data_entities
    assert "hasPart" not in crate.root_dataset


def test_entity_as_mapping(tmpdir, helpers):
    orcid = "https://orcid.org/0000-0002-1825-0097"
    metadata = {
        "@context": "https://w3id.org/ro/crate/1.1/context",
        "@graph": [
            {"@id": "ro-crate-metadata.json",
             "@type": "CreativeWork",
             "about": {"@id": "./"},
             "encodingFormat": [
                 "application/json",
                 {"@id": "https://www.json.org"},
             ],
             "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"}},
            {"@id": "./",
             "@type": "Dataset",
             "correction": [
                 "Fixed typo.",
                 {"@id": "#correction"},
                 "http://example.org/correction",
             ],
             "author": {"@id": orcid}},
            {"@id": "#correction",
             "@type": "CorrectionComment",
             "badProp": {"k": "v"},
             "text": "Previous version was ugly."},
            {"@id": orcid,
             "@type": "Person",
             "name": None,
             "givenName": "Josiah",
             "familyName": "Carberry"}
        ]
    }
    crate_dir = tmpdir / "in_crate"
    crate_dir.mkdir()
    with open(crate_dir / helpers.METADATA_FILE_NAME, "wt") as f:
        json.dump(metadata, f, indent=4)
    with pytest.raises(ValueError):  # due to "badProp", which has no "@id"
        crate = ROCrate(crate_dir)
    del metadata["@graph"][2]["badProp"]
    with open(crate_dir / helpers.METADATA_FILE_NAME, "wt") as f:
        json.dump(metadata, f, indent=4)
    crate = ROCrate(crate_dir)
    person = crate.dereference(orcid)
    exp_len = len([_ for _ in metadata["@graph"] if _["@id"] == orcid][0])
    assert len(person) == exp_len
    assert len(list(person)) == exp_len
    assert set(person) == set(person.keys()) == {"@id", "@type", "name", "givenName", "familyName"}
    assert set(person.values()) == {orcid, "Person", None, "Josiah", "Carberry"}
    assert set(person.items()) == set(zip(person.keys(), person.values()))
    assert person.id == orcid
    assert person.type == "Person"
    assert person["name"] is None
    assert person["givenName"] == "Josiah"
    assert person["familyName"] == "Carberry"
    assert person.setdefault("givenName", "foo") == "Josiah"
    assert len(person) == exp_len
    assert person.setdefault("award", "Oscar") == "Oscar"
    assert len(person) == exp_len + 1
    assert "award" in person
    assert person.pop("award") == "Oscar"
    assert len(person) == exp_len
    assert "award" not in person
    for key in "@id", "@type":
        with pytest.raises(KeyError):
            person[key] = "foo"
        with pytest.raises(KeyError):
            del person[key]
        with pytest.raises(KeyError):
            person.pop(key)
    twin = Person(crate, orcid, properties={
        "name": None,
        "givenName": "Josiah",
        "familyName": "Carberry"
    })
    assert twin == person
    assert Person(crate, orcid) != person
    assert crate.root_dataset["author"] is person
    correction = crate.get("#correction")
    assert set(crate.root_dataset["correction"]) == {
        "Fixed typo.",
        correction,
        "http://example.org/correction"
    }
    assert set(crate.metadata["encodingFormat"]) == {
        "application/json",
        "https://www.json.org",
    }
    with pytest.raises(ValueError):
        correction["badProp"] = {"k": "v"}  # value has no "@id"
    correction._jsonld["badProp"] = {"k": "v"}  # force set using _jsonld
    with pytest.raises(ValueError):
        correction["badProp"]


def test_wf_types():
    foo_crate = ROCrate()
    foo_wf = foo_crate.add_workflow("foo.cwl", main=True)
    assert "HowTo" not in foo_wf.type
    foo_wf.type.append("HowTo")
    bar_crate = ROCrate()
    bar_wf = bar_crate.add_workflow("bar.cwl", main=True)
    assert "HowTo" not in bar_wf.type


@pytest.mark.parametrize("compact", [False, True])
def test_append_to(compact):
    crate = ROCrate()
    alice = crate.add(Person(crate, "#alice"))
    bob = crate.add(Person(crate, "#bob"))
    rd = crate.root_dataset
    assert rd.get("author") is None
    rd.append_to("author", alice, compact=compact)
    assert rd.get("author") == (alice if compact else [alice])
    rd.append_to("author", bob, compact=compact)
    assert set(rd.get("author")) == {alice, bob}
    # as string
    don = "https://en.wikipedia.org/wiki/Donald_Duck"
    rd.append_to("author", don, compact=compact)
    assert set(rd.get("author")) == {alice, bob, don}
    # multiple values
    scrooge = "https://en.wikipedia.org/wiki/Scrooge_McDuck"
    charlie = crate.add(Person(crate, "#charlie"))
    rd.append_to("author", [scrooge, charlie], compact=compact)
    assert set(rd.get("author")) == {alice, bob, don, scrooge, charlie}
    # exceptions
    with pytest.raises(KeyError):
        rd.append_to("@id", "foo", compact=compact)


def test_get_by_type(test_data_dir):
    crate = ROCrate()
    f1 = crate.add_file(test_data_dir / "sample_file.txt")
    f2 = crate.add_file(test_data_dir / "test_file_galaxy.txt")
    wf = crate.add_workflow(test_data_dir / "test_galaxy_wf.ga", main=True, lang="galaxy")
    assert set(crate.get_by_type("File")) == {f1, f2, wf}
    assert set(crate.get_by_type("File", exact=True)) == {f1, f2}
    assert crate.get_by_type(["ComputationalWorkflow"]) == [wf]
    assert crate.get_by_type(["File", "ComputationalWorkflow"]) == [wf]
    assert crate.get_by_type("Person") == []
    assert crate.get_by_type("Dataset") == [crate.root_dataset]


def test_context(helpers):
    crate = ROCrate()
    jsonld = crate.metadata.generate()
    base_context = f"{helpers.PROFILE}/context"
    assert jsonld["@context"] == base_context
    wfrun_ctx = "https://w3id.org/ro/terms/workflow-run"
    crate.metadata.extra_contexts.append(wfrun_ctx)
    jsonld = crate.metadata.generate()
    assert jsonld["@context"] == [base_context, wfrun_ctx]
    k, v = "runsOn", "https://w3id.org/ro/terms/test#runsOn"
    crate.metadata.extra_terms[k] = v
    jsonld = crate.metadata.generate()
    assert jsonld["@context"] == [base_context, wfrun_ctx, {k: v}]


def test_add_no_duplicates(test_data_dir, tmpdir):
    source = test_data_dir / "sample_file.txt"
    crate = ROCrate()
    f1 = crate.add_file(source, properties={"name": "sample file"})
    ret = crate.get(source.name)
    assert ret is f1
    assert ret["name"] == "sample file"
    assert ret in crate.get_entities()
    assert crate.data_entities == [f1]
    f2 = crate.add_file(source, properties={"name": "foobar"})
    ret = crate.get(source.name)
    assert ret is f2
    assert ret["name"] == "foobar"
    assert ret in crate.get_entities()
    assert f1 not in crate.get_entities()
    assert crate.data_entities == [f2]
    joe = crate.add(Person(crate, "#joe", properties={"name": "Joe"}))
    ret = crate.get("#joe")
    assert ret is joe
    assert ret in crate.get_entities()
    assert ret["name"] == "Joe"
    assert crate.contextual_entities == [joe]
    jim = crate.add(Person(crate, "#joe", properties={"name": "Jim"}))
    ret = crate.get("#joe")
    assert ret is jim
    assert ret["name"] == "Jim"
    assert ret in crate.get_entities()
    assert crate.contextual_entities == [jim]


def test_immutable_id():
    crate = ROCrate()
    p = crate.add(Person(crate, "#foo"))
    with pytest.raises(AttributeError):
        p.id = "#bar"


def test_entity_in_properties(tmpdir):
    file_id = "data.csv"
    file_path = tmpdir / file_id
    with open(file_path, "wt") as fout:
        fout.write("DATA\n")
    crate = ROCrate()
    alice_id = "https://orcid.org/0000-0000-0000-0000"
    bob_id = "https://orcid.org/0000-0000-0000-0001"
    alice = crate.add(Person(crate, alice_id, properties={
        "name": "Alice Doe",
        "affiliation": "University of Flatland",
    }))
    bob = crate.add(Person(crate, bob_id, properties={
        "name": "Bob Doe",
        "affiliation": "University of Flatland",
    }))
    data = crate.add_file(file_path, properties={
        "name": "Data file",
        "encodingFormat": "text/csv",
        "author": [alice, bob],
    })
    authors = data["author"]
    assert len(authors) == 2
    assert authors[0] is alice
    assert authors[1] is bob

    book = crate.add(ContextEntity(crate, "#44cats", properties={
        "@type": "Book",
        "name": "44 Cats",
        "author": alice,
    }))
    assert book.type == "Book"
    assert book["author"] is alice

    out_path = tmpdir / "ro_crate_out"
    crate.write(out_path)
    assert (out_path / file_id).is_file()
    out_crate = ROCrate(out_path)
    out_alice = out_crate.get(alice_id)
    out_bob = out_crate.get(bob_id)
    out_data = out_crate.get(file_id)
    out_authors = out_data["author"]
    assert len(out_authors) == 2
    assert out_authors[0] is out_alice
    assert out_authors[1] is out_bob
    assert out_alice == alice
    assert out_bob == bob
