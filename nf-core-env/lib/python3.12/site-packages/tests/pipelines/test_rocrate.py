"""Test the nf-core pipelines rocrate command"""

import json
import shutil
import tempfile
from pathlib import Path

import git
import rocrate.rocrate
from git import Repo

import nf_core.pipelines.create
import nf_core.pipelines.create.create
import nf_core.pipelines.rocrate
import nf_core.utils
from nf_core.pipelines.bump_version import bump_pipeline_version

from ..test_pipelines import TestPipelines


class TestROCrate(TestPipelines):
    """Class for lint tests"""

    def setUp(self) -> None:
        super().setUp()
        # add fake metro map
        Path(self.pipeline_dir, "docs", "images", "nf-core-testpipeline_metro_map.png").touch()
        # commit the changes
        repo = Repo(self.pipeline_dir)
        repo.git.add(A=True)
        repo.index.commit("Initial commit")
        self.rocrate_obj = nf_core.pipelines.rocrate.ROCrate(self.pipeline_dir)

    def tearDown(self):
        """Clean up temporary files and folders"""

        if self.tmp_dir.exists():
            shutil.rmtree(self.tmp_dir)

    def test_rocrate_creation(self):
        """Run the nf-core rocrate command"""

        # Run the command
        self.rocrate_obj
        assert self.rocrate_obj.create_rocrate(self.pipeline_dir, self.pipeline_dir)

        # Check that the crate was created
        self.assertTrue(Path(self.pipeline_dir, "ro-crate-metadata.json").exists())

        # Check that the entries in the crate are correct
        crate = rocrate.rocrate.ROCrate(self.pipeline_dir)
        entities = crate.get_entities()

        # Check if the correct entities are set:
        for entity in entities:
            entity_json = entity.as_jsonld()
            if entity_json["@id"] == "./":
                self.assertEqual(
                    entity_json.get("name"), f"{self.pipeline_obj.pipeline_prefix}/{self.pipeline_obj.pipeline_name}"
                )
                self.assertEqual(entity_json["mainEntity"], {"@id": "main.nf"})
            elif entity_json["@id"] == "#main.nf":
                self.assertEqual(entity_json["programmingLanguage"], [{"@id": "#nextflow"}])
                self.assertEqual(entity_json["image"], [{"@id": "nf-core-testpipeline_metro_map.png"}])
            # assert there is a metro map
            # elif entity_json["@id"] == "nf-core-testpipeline_metro_map.png": # FIXME waiting for https://github.com/ResearchObject/ro-crate-py/issues/174
            # self.assertEqual(entity_json["@type"], ["File", "ImageObject"])
            # assert that author is set as a person
            elif "name" in entity_json and entity_json["name"] == "Test McTestFace":
                self.assertEqual(entity_json["@type"], "Person")
                # check that it is set as author of the main entity
                if crate.mainEntity is not None:
                    self.assertEqual(crate.mainEntity["author"][0].id, entity_json["@id"])

    def test_rocrate_creation_wrong_pipeline_dir(self):
        """Run the nf-core rocrate command with a wrong pipeline directory"""
        # Run the command

        # Check that it raises a UserWarning
        with self.assertRaises(UserWarning):
            nf_core.pipelines.rocrate.ROCrate(self.pipeline_dir / "bad_dir")

        # assert that the crate was not created
        self.assertFalse(Path(self.pipeline_dir / "bad_dir", "ro-crate-metadata.json").exists())

    def test_rocrate_creation_with_wrong_version(self):
        """Run the nf-core rocrate command with a pipeline version"""
        # Run the command

        self.rocrate_obj = nf_core.pipelines.rocrate.ROCrate(self.pipeline_dir, version="1.0.0")

        # Check that the crate was created
        with self.assertRaises(SystemExit):
            assert self.rocrate_obj.create_rocrate(self.pipeline_dir, self.pipeline_dir)

    def test_rocrate_creation_without_git(self):
        """Run the nf-core rocrate command with a pipeline version"""

        self.rocrate_obj = nf_core.pipelines.rocrate.ROCrate(self.pipeline_dir, version="1.0.0")
        # remove git repo
        shutil.rmtree(self.pipeline_dir / ".git")
        # Check that the crate was created
        with self.assertRaises(SystemExit):
            assert self.rocrate_obj.create_rocrate(self.pipeline_dir, self.pipeline_dir)

    def test_rocrate_creation_to_zip(self):
        """Run the nf-core rocrate command with a zip output"""
        assert self.rocrate_obj.create_rocrate(self.pipeline_dir, zip_path=self.pipeline_dir)
        # Check that the crate was created
        self.assertTrue(Path(self.pipeline_dir, "ro-crate.crate.zip").exists())

    def test_rocrate_creation_for_fetchngs(self):
        """Run the nf-core rocrate command with nf-core/fetchngs"""
        tmp_dir = Path(tempfile.mkdtemp())
        # git clone  nf-core/fetchngs
        git.Repo.clone_from("https://github.com/nf-core/fetchngs", tmp_dir / "fetchngs")
        # Run the command
        self.rocrate_obj = nf_core.pipelines.rocrate.ROCrate(tmp_dir / "fetchngs", version="1.12.0")
        assert self.rocrate_obj.create_rocrate(tmp_dir / "fetchngs", self.pipeline_dir)

        # Check that Sateesh Peri is mentioned in creator field

        crate = rocrate.rocrate.ROCrate(self.pipeline_dir)
        entities = crate.get_entities()
        for entity in entities:
            entity_json = entity.as_jsonld()
            if entity_json["@id"] == "#main.nf":
                assert "https://orcid.org/0000-0002-9879-9070" in entity_json["creator"]

        # Clean up
        shutil.rmtree(tmp_dir)

    def test_update_rocrate(self):
        """Run the nf-core rocrate command with a zip output"""

        assert self.rocrate_obj.create_rocrate(json_path=self.pipeline_dir, zip_path=self.pipeline_dir)

        # read the crate json file
        with open(Path(self.pipeline_dir, "ro-crate-metadata.json")) as f:
            crate = json.load(f)

        # check the old version
        self.assertEqual(crate["@graph"][2]["version"][0], "1.0.0dev")
        # check creativeWorkStatus is InProgress
        self.assertEqual(crate["@graph"][0]["creativeWorkStatus"], "InProgress")

        # bump version
        bump_pipeline_version(self.pipeline_obj, "1.1.0")

        # Check that the crate was created
        self.assertTrue(Path(self.pipeline_dir, "ro-crate.crate.zip").exists())

        # Check that the crate was updated
        self.assertTrue(Path(self.pipeline_dir, "ro-crate-metadata.json").exists())

        # read the crate json file
        with open(Path(self.pipeline_dir, "ro-crate-metadata.json")) as f:
            crate = json.load(f)

        # check that the version was updated
        self.assertEqual(crate["@graph"][2]["version"][0], "1.1.0")

        # check creativeWorkStatus is Stable
        self.assertEqual(crate["@graph"][0]["creativeWorkStatus"], "Stable")
