#!/usr/bin/env python
"""Code to deal with pipeline RO (Research Object) Crates"""

import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import requests
import rocrate.rocrate
from git import GitCommandError, InvalidGitRepositoryError
from repo2rocrate.nextflow import NextflowCrateBuilder
from rich.progress import BarColumn, Progress
from rocrate.model.person import Person
from rocrate.rocrate import ROCrate as BaseROCrate

from nf_core.utils import Pipeline

log = logging.getLogger(__name__)


class CustomNextflowCrateBuilder(NextflowCrateBuilder):
    DATA_ENTITIES = NextflowCrateBuilder.DATA_ENTITIES + [
        ("docs/usage.md", "File", "Usage documentation"),
        ("docs/output.md", "File", "Output documentation"),
        ("suborkflows/local", "Dataset", "Pipeline-specific suborkflows"),
        ("suborkflows/nf-core", "Dataset", "nf-core suborkflows"),
        (".nf-core.yml", "File", "nf-core configuration file, configuring template features and linting rules"),
        (".pre-commit-config.yaml", "File", "Configuration file for pre-commit hooks"),
        (".prettierignore", "File", "Ignore file for prettier"),
        (".prettierrc", "File", "Configuration file for prettier"),
    ]


def custom_make_crate(
    root: Path,
    workflow: Path | None = None,
    repo_url: str | None = None,
    wf_name: str | None = None,
    wf_version: str | None = None,
    lang_version: str | None = None,
    ci_workflow: str | None = "nf-test.yml",
    diagram: Path | None = None,
) -> BaseROCrate:
    builder = CustomNextflowCrateBuilder(root, repo_url=repo_url)

    return builder.build(
        workflow,
        wf_name=wf_name,
        wf_version=wf_version,
        lang_version=lang_version,
        license=None,
        ci_workflow=ci_workflow,
        diagram=diagram,
    )


class ROCrate:
    """
    Class to generate an RO Crate for a pipeline

    """

    def __init__(self, pipeline_dir: Path, version="") -> None:
        """
        Initialise the ROCrate object

        Args:
            pipeline_dir (Path): Path to the pipeline directory
            version (str): Version of the pipeline to checkout
        """
        from nf_core.utils import is_pipeline_directory, setup_requests_cachedir

        is_pipeline_directory(pipeline_dir)
        self.pipeline_dir = pipeline_dir
        self.version: str = version
        self.crate: rocrate.rocrate.ROCrate
        self.pipeline_obj = Pipeline(self.pipeline_dir)
        self.pipeline_obj._load()

        setup_requests_cachedir()

    def create_rocrate(self, json_path: None | Path = None, zip_path: None | Path = None) -> bool:
        """
        Create an RO Crate for a pipeline

        Args:
            outdir (Path): Path to the output directory
            json_path (Path): Path to the metadata file
            zip_path (Path): Path to the zip file

        """

        # Check that the checkout pipeline version is the same as the requested version
        if self.version != "":
            if self.version != self.pipeline_obj.nf_config.get("manifest.version"):
                # using git checkout to get the requested version
                log.info(f"Checking out pipeline version {self.version}")
                if self.pipeline_obj.repo is None:
                    log.error(f"Pipeline repository not found in {self.pipeline_dir}")
                    sys.exit(1)
                try:
                    self.pipeline_obj.repo.git.checkout(self.version)
                    self.pipeline_obj = Pipeline(self.pipeline_dir)
                    self.pipeline_obj._load()
                except InvalidGitRepositoryError:
                    log.error(f"Could not find a git repository in {self.pipeline_dir}")
                    sys.exit(1)
                except GitCommandError:
                    log.error(f"Could not checkout version {self.version}")
                    sys.exit(1)
        self.version = self.pipeline_obj.nf_config.get("manifest.version", "")
        self.make_workflow_rocrate()

        # Save just the JSON metadata file
        if json_path is not None:
            if json_path.name == "ro-crate-metadata.json":
                json_path = json_path.parent

            log.info(f"Saving metadata file to '{json_path}'")
            self.crate.metadata.write(json_path)

        # Save the whole crate zip file
        if zip_path is not None:
            if zip_path.name != "ro-crate.crate.zip":
                zip_path = zip_path / "ro-crate.crate.zip"
            log.info(f"Saving zip file '{zip_path}")
            self.crate.write_zip(zip_path)

        if json_path is None and zip_path is None:
            log.error("Please provide a path to save the ro-crate file or the zip file.")
            return False

        return True

    def make_workflow_rocrate(self) -> None:
        """
        Create an RO Crate for a pipeline
        """
        if self.pipeline_obj is None:
            raise ValueError("Pipeline object not loaded")

        diagram: Path | None = None
        # find files (metro|tube)_?(map)?.png in the pipeline directory or docs/ using pathlib
        pattern = re.compile(r".*?(metro|tube|subway)_(map).*?\.png", re.IGNORECASE)
        for file in self.pipeline_dir.rglob("*.png"):
            if pattern.match(file.name):
                log.debug(f"Found diagram: {file}")
                diagram = file.relative_to(self.pipeline_dir)
                break

        # Create the RO Crate object

        self.crate = custom_make_crate(
            self.pipeline_dir,
            self.pipeline_dir / "main.nf",
            self.pipeline_obj.nf_config.get("manifest.homePage", ""),
            self.pipeline_obj.nf_config.get("manifest.name", ""),
            self.pipeline_obj.nf_config.get("manifest.version", ""),
            self.pipeline_obj.nf_config.get("manifest.nextflowVersion", ""),
            diagram=diagram,
        )

        # add readme as description
        readme = self.pipeline_dir / "README.md"

        try:
            self.crate.description = readme.read_text()
        except FileNotFoundError:
            log.error(f"Could not find README.md in {self.pipeline_dir}")
        # get license from LICENSE file
        license_file = self.pipeline_dir / "LICENSE"
        try:
            license = license_file.read_text()
            if license.startswith("MIT"):
                self.crate.license = "MIT"
            else:
                # prompt for license
                log.info("Could not determine license from LICENSE file")
                self.crate.license = input("Please enter the license for this pipeline: ")
        except FileNotFoundError:
            log.error(f"Could not find LICENSE file in {self.pipeline_dir}")

        self.crate.add_jsonld(
            {"@id": "https://nf-co.re/", "@type": "Organization", "name": "nf-core", "url": "https://nf-co.re/"}
        )

        # Set metadata for main entity file
        self.set_main_entity("main.nf")

    def set_main_entity(self, main_entity_filename: str):
        """
        Set the main.nf as the main entity of the crate and add necessary metadata
        """
        if self.crate.mainEntity is None:
            raise ValueError("Main entity not set")

        self.crate.mainEntity.append_to(
            "dct:conformsTo", "https://bioschemas.org/profiles/ComputationalWorkflow/1.0-RELEASE/", compact=True
        )
        # add dateCreated and dateModified, based on the current data
        self.crate.mainEntity.append_to("dateCreated", self.crate.root_dataset.get("dateCreated", ""), compact=True)
        self.crate.mainEntity.append_to(
            "dateModified", str(datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")), compact=True
        )
        self.crate.mainEntity.append_to("sdPublisher", {"@id": "https://nf-co.re/"}, compact=True)
        if self.version.endswith("dev"):
            url = "dev"
        else:
            url = self.version
        self.crate.mainEntity.append_to(
            "url", f"https://nf-co.re/{self.crate.name.replace('nf-core/', '')}/{url}/", compact=True
        )
        self.crate.mainEntity.append_to("version", self.version, compact=True)

        # remove duplicate entries for version
        self.crate.mainEntity["version"] = list(set(self.crate.mainEntity["version"]))

        # get keywords from nf-core website
        remote_workflows = requests.get("https://nf-co.re/pipelines.json").json()["remote_workflows"]
        # go through all remote workflows and find the one that matches the pipeline name
        topics = ["nf-core", "nextflow"]
        for remote_wf in remote_workflows:
            assert self.pipeline_obj.pipeline_name is not None  # mypy
            if remote_wf["name"] == self.pipeline_obj.pipeline_name.replace("nf-core/", ""):
                topics = topics + remote_wf["topics"]
                break

        log.debug(f"Adding topics: {topics}")
        self.crate.mainEntity.append_to("keywords", topics)

        self.add_main_authors(self.crate.mainEntity)

        self.crate.mainEntity = self.crate.mainEntity

        self.crate.mainEntity.append_to("license", self.crate.license)
        self.crate.mainEntity.append_to("name", self.crate.name)

        # remove duplicate entries for name
        self.crate.mainEntity["name"] = list(set(self.crate.mainEntity["name"]))

        if "dev" in self.version:
            self.crate.creativeWorkStatus = "InProgress"
        else:
            self.crate.creativeWorkStatus = "Stable"
            if self.pipeline_obj.repo is None:
                log.error(f"Pipeline repository not found in {self.pipeline_dir}")
            else:
                tags = self.pipeline_obj.repo.tags
                if tags:
                    # get the tag for this version
                    for tag in tags:
                        if tag.commit.hexsha == self.pipeline_obj.repo.head.commit.hexsha:
                            self.crate.mainEntity.append_to(
                                "dateCreated",
                                tag.commit.committed_datetime.strftime("%Y-%m-%dT%H:%M:%SZ"),
                                compact=True,
                            )

    def add_main_authors(self, wf_file: rocrate.model.entity.Entity) -> None:
        """
        Add workflow authors to the crate
        """
        # add author entity to crate

        try:
            authors = []
            if "manifest.author" in self.pipeline_obj.nf_config:
                authors.extend([a.strip() for a in self.pipeline_obj.nf_config["manifest.author"].split(",")])
            if "manifest.contributors" in self.pipeline_obj.nf_config:
                contributors = self.pipeline_obj.nf_config["manifest.contributors"]
                names = re.findall(r"name:'([^']+)'", contributors)
                authors.extend(names)
            if not authors:
                raise KeyError("No authors found")
            # add manifest authors as maintainer to crate

        except KeyError:
            log.error("No author or contributors fields found in manifest of nextflow.config")
            return
        # remove duplicates
        authors = list(set(authors))
        # look at git contributors for author names
        try:
            git_contributors: set[str] = set()
            if self.pipeline_obj.repo is None:
                log.debug("No git repository found. No git contributors will be added as authors.")
                return
            commits_touching_path = list(self.pipeline_obj.repo.iter_commits(paths="main.nf"))

            for commit in commits_touching_path:
                if commit.author.name is not None:
                    git_contributors.add(commit.author.name)
            # exclude bots
            contributors = {c for c in git_contributors if not c.endswith("bot") and c != "Travis CI User"}

            log.debug(f"Found {len(contributors)} git authors")

            progress_bar = Progress(
                "[bold blue]{task.description}",
                BarColumn(bar_width=None),
                "[magenta]{task.completed} of {task.total}[reset] » [bold yellow]{task.fields[test_name]}",
                transient=True,
                disable=os.environ.get("HIDE_PROGRESS", None) is not None,
            )
            with progress_bar:
                bump_progress = progress_bar.add_task(
                    "Searching for author names on GitHub", total=len(contributors), test_name=""
                )

                for git_author in contributors:
                    progress_bar.update(bump_progress, advance=1, test_name=git_author)
                    git_author = (
                        requests.get(f"https://api.github.com/users/{git_author}").json().get("name", git_author)
                    )
                    if git_author is None:
                        log.debug(f"Could not find name for {git_author}")
                        continue

        except AttributeError:
            log.debug("Could not find git contributors")

        # remove usernames (just keep names with spaces)
        named_contributors = {c for c in contributors if " " in c}

        for author in named_contributors:
            log.debug(f"Adding author: {author}")

            if self.pipeline_obj.repo is None:
                log.info("No git repository found. No git contributors will be added as authors.")
                return
            # get email from git log
            email = self.pipeline_obj.repo.git.log(f"--author={author}", "--pretty=format:%ae", "-1")
            orcid = get_orcid(author)
            author_entitity = self.crate.add(
                Person(
                    self.crate, orcid if orcid is not None else "#" + email, properties={"name": author, "email": email}
                )
            )
            wf_file.append_to("creator", author_entitity)
            if author in authors:
                wf_file.append_to("maintainer", author_entitity)

    def update_rocrate(self) -> bool:
        """
        Update the rocrate file
        """
        # check if we need to output a json file and/or a zip file based on the file extensions
        # try to find a json file
        json_path: Path | None = None
        potential_json_path = Path(self.pipeline_dir, "ro-crate-metadata.json")
        if potential_json_path.exists():
            json_path = potential_json_path

        # try to find a zip file
        zip_path: Path | None = None
        potential_zip_path = Path(self.pipeline_dir, "ro-crate.crate.zip")
        if potential_zip_path.exists():
            zip_path = potential_zip_path

        return self.create_rocrate(json_path=json_path, zip_path=zip_path)


def get_orcid(name: str) -> str | None:
    """
    Get the ORCID for a given name

    Args:
        name (str): Name of the author

    Returns:
        str: ORCID URI or None
    """
    base_url = "https://pub.orcid.org/v3.0/search/"
    headers = {
        "Accept": "application/json",
    }
    params = {"q": f'family-name:"{name.split()[-1]}" AND given-names:"{name.split()[0]}"'}
    response = requests.get(base_url, params=params, headers=headers)

    if response.status_code == 200:
        json_response = response.json()
        if json_response.get("num-found") == 1:
            orcid_uri = json_response.get("result")[0].get("orcid-identifier", {}).get("uri")
            log.info(f"Using found ORCID for {name}. Please double-check: {orcid_uri}")
            return orcid_uri
        else:
            log.debug(f"No exact ORCID found for {name}. See {response.url}")
            return None
    else:
        log.info(f"API request to ORCID unsuccessful. Status code: {response.status_code}")
        return None
