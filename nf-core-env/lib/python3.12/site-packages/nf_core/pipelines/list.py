"""Lists available nf-core pipelines and versions."""

import json
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import git
import requests
import rich.console
import rich.table
from click.shell_completion import CompletionItem

import nf_core.utils

log = logging.getLogger(__name__)

# Set up local caching for requests to speed up remote queries
nf_core.utils.setup_requests_cachedir()


def list_workflows(filter_by=None, sort_by="release", as_json=False, show_archived=False):
    """Prints out a list of all nf-core workflows.

    Args:
        filter_by (list): A list of strings that can be used for filtering.
        sort_by (str): workflows can be sorted by keywords. Keyword must be one of
            `release` (default), `name`, `stars`.
        as_json (boolean): Set to true, if the lists should be printed in JSON.
    """
    wfs = Workflows(filter_by, sort_by, show_archived)
    wfs.get_remote_workflows()
    wfs.get_local_nf_workflows()
    wfs.compare_remote_local()
    if as_json:
        return wfs.print_json()
    else:
        return wfs.print_summary()


def autocomplete_pipelines(ctx, param, incomplete: str):
    try:
        wfs = Workflows()
        wfs.get_remote_workflows()
        wfs.get_local_nf_workflows()
        local_workflows = [wf.full_name for wf in wfs.local_workflows]
        remote_workflows = [wf.full_name for wf in wfs.remote_workflows]
        available_workflows = local_workflows + remote_workflows

        matches = [CompletionItem(wor) for wor in available_workflows if wor.startswith(incomplete)]

        return matches
    except Exception as e:
        print(f"[ERROR] Autocomplete failed: {e}", file=sys.stderr)
        return []


def get_local_wf(workflow: str | Path, revision=None) -> str | None:
    """
    Check if this workflow has a local copy and use nextflow to pull it if not
    """
    # Assume nf-core if no org given
    if str(workflow).count("/") == 0:
        workflow = f"nf-core/{workflow}"

    wfs = Workflows()
    wfs.get_local_nf_workflows()
    for wf in wfs.local_workflows:
        if workflow == wf.full_name:
            if revision is None or revision == wf.commit_sha or revision == wf.branch or revision == wf.active_tag:
                if wf.active_tag:
                    print_revision = f"v{wf.active_tag}"
                elif wf.branch:
                    print_revision = f"{wf.branch} - {wf.commit_sha[:7]}"
                else:
                    print_revision = wf.commit_sha
                log.info(f"Using local workflow: {workflow} ({print_revision})")
                return wf.local_path

    # Wasn't local, fetch it
    log.info(f"Downloading workflow: {workflow} ({revision})")
    pull_cmd = f"pull {workflow}"
    if revision is not None:
        pull_cmd += f" -r {revision}"
    nf_core.utils.run_cmd("nextflow", pull_cmd)
    local_wf = LocalWorkflow(workflow)
    local_wf.get_local_nf_workflow_details()
    return local_wf.local_path


class Workflows:
    """Workflow container class.

    Is used to collect local and remote nf-core pipelines. Pipelines
    can be sorted, filtered and compared.

    Args:
        filter_by (list): A list of strings that can be used for filtering.
        sort_by (str): workflows can be sorted by keywords. Keyword must be one of
            `release` (default), `name`, `stars`.
    """

    def __init__(self, filter_by=None, sort_by="release", show_archived=False):
        self.remote_workflows = []
        self.local_workflows = []
        self.local_unmatched = []
        self.keyword_filters = filter_by if filter_by is not None else []
        self.sort_workflows_by = sort_by
        self.show_archived = show_archived

    def get_remote_workflows(self):
        """Retrieves remote workflows from `nf-co.re <https://nf-co.re>`_.

        Remote workflows are stored in :attr:`self.remote_workflows` list.
        """
        # List all repositories at nf-core
        log.debug("Fetching list of nf-core workflows")
        nfcore_url = "https://nf-co.re/pipelines.json"
        response = requests.get(nfcore_url, timeout=10)
        if response.status_code == 200:
            repos = response.json()["remote_workflows"]
            for repo in repos:
                self.remote_workflows.append(RemoteWorkflow(repo))

    def get_local_nf_workflows(self):
        """Retrieves local Nextflow workflows.

        Local workflows are stored in :attr:`self.local_workflows` list.
        """
        # Try to guess the local cache directory (much faster than calling nextflow)
        if len(os.environ.get("NXF_ASSETS", "")) > 0:
            nextflow_wfdir = os.environ.get("NXF_ASSETS")
        elif len(os.environ.get("NXF_HOME", "")) > 0:
            nextflow_wfdir = os.path.join(os.environ.get("NXF_HOME"), "assets")
        else:
            nextflow_wfdir = os.path.join(os.getenv("HOME"), ".nextflow", "assets")
        if os.path.isdir(nextflow_wfdir):
            log.debug("Guessed nextflow assets directory - pulling pipeline dirnames")
            for org_name in os.listdir(nextflow_wfdir):
                for wf_name in os.listdir(os.path.join(nextflow_wfdir, org_name)):
                    self.local_workflows.append(LocalWorkflow(f"{org_name}/{wf_name}"))

        # Fetch details about local cached pipelines with `nextflow list`
        else:
            log.debug("Getting list of local nextflow workflows")
            result = nf_core.utils.run_cmd("nextflow", "list")
            if result is not None:
                nflist_raw, _ = result
                for wf_name in nflist_raw.splitlines():
                    if not str(wf_name).startswith("nf-core/"):
                        self.local_unmatched.append(wf_name)
                    else:
                        self.local_workflows.append(LocalWorkflow(wf_name))

        # Find additional information about each workflow by checking its git history
        log.debug(f"Fetching extra info about {len(self.local_workflows)} local workflows")
        for wf in self.local_workflows:
            wf.get_local_nf_workflow_details()

    def compare_remote_local(self):
        """Matches local to remote workflows.

        If a matching remote workflow is found, the local workflow's Git commit hash is compared
        with the latest one from remote.

        A boolean flag in :attr:`RemoteWorkflow.local_is_latest` is set to True, if the local workflow
        is the latest.
        """
        for rwf in self.remote_workflows:
            for lwf in self.local_workflows:
                if rwf.full_name == lwf.full_name:
                    rwf.local_wf = lwf
                    if rwf.releases:
                        if rwf.releases[-1]["tag_sha"] == lwf.commit_sha:
                            rwf.local_is_latest = True
                        else:
                            rwf.local_is_latest = False

    def filtered_workflows(self):
        """Filters remote workflows for keywords.

        Returns:
            list: Filtered remote workflows.
        """
        filtered_workflows = []
        for wf in self.remote_workflows:
            # Skip archived pipelines
            if not self.show_archived and wf.archived:
                continue
            # Search through any supplied keywords
            for k in self.keyword_filters:
                in_name = k in wf.name if wf.name else False
                in_desc = k in wf.description if wf.description else False
                in_topics = any(k in t for t in wf.topics)
                if not in_name and not in_desc and not in_topics:
                    break
            else:
                # We didn't hit a break, so all keywords were found
                filtered_workflows.append(wf)

        return filtered_workflows

    def print_summary(self):
        """Prints a summary of all pipelines."""

        filtered_workflows = self.filtered_workflows()

        # Sort by released / dev / archived, then alphabetical
        if not self.sort_workflows_by or self.sort_workflows_by == "release":
            filtered_workflows.sort(
                key=lambda wf: (
                    (wf.releases[-1].get("published_at_timestamp", 0) if len(wf.releases) > 0 else 0) * -1,
                    wf.archived,
                    wf.full_name.lower(),
                )
            )
        # Sort by date pulled
        elif self.sort_workflows_by == "pulled":

            def sort_pulled_date(wf):
                try:
                    return wf.local_wf.last_pull * -1
                except Exception:
                    return 0

            filtered_workflows.sort(key=sort_pulled_date)
        # Sort by name
        elif self.sort_workflows_by == "name":
            filtered_workflows.sort(key=lambda wf: wf.full_name.lower())
        # Sort by stars, then name
        elif self.sort_workflows_by == "stars":
            filtered_workflows.sort(key=lambda wf: (wf.stargazers_count * -1, wf.full_name.lower()))

        # Build summary list to print
        table = rich.table.Table()
        table.add_column("Pipeline Name")
        table.add_column("Stars", justify="right")
        table.add_column("Latest Release", justify="right")
        table.add_column("Released", justify="right")
        table.add_column("Last Pulled", justify="right")
        table.add_column("Have latest release?")
        for wf in filtered_workflows:
            wf_name = f"[bold][link=https://nf-co.re/{wf.name}]{wf.name}[/link]"
            version = "[yellow]dev"
            published = "[dim]-"
            if len(wf.releases) > 0:
                version = f"[blue]{wf.releases[0]['tag_name']}"
                published = wf.releases[0]["published_at_pretty"]
            pulled = wf.local_wf.last_pull_pretty if wf.local_wf is not None else "[dim]-"
            if wf.local_wf is not None:
                revision = ""
                if wf.local_wf.active_tag is not None:
                    revision = f"v{wf.local_wf.active_tag}"
                elif wf.local_wf.branch is not None:
                    revision = f"{wf.local_wf.branch} - {wf.local_wf.commit_sha[:7]}"
                else:
                    revision = wf.local_wf.commit_sha
                if wf.local_is_latest:
                    is_latest = f"[green]Yes ({revision})"
                else:
                    is_latest = f"[red]No ({revision})"
            else:
                is_latest = "[dim]-"

            rowdata = [wf_name, str(wf.stargazers_count), version, published, pulled, is_latest]

            # Handle archived pipelines
            if wf.archived:
                rowdata[1] = "archived"
                rowdata = [re.sub(r"\[\w+\]", "", k) for k in rowdata]
                table.add_row(*rowdata, style="dim")
            else:
                table.add_row(*rowdata)

        if len(filtered_workflows) > 0:
            # Print summary table
            return table
        else:
            return_str = f"No pipelines found using filter keywords: '{', '.join(self.keyword_filters)}'"
            if self.keyword_filters == ("modules",):
                return_str += "\n\n:bulb: Did you mean 'nf-core modules list' instead?"
            return return_str

    def print_json(self):
        """Dump JSON of all parsed information"""
        return json.dumps(
            {"local_workflows": self.local_workflows, "remote_workflows": self.remote_workflows},
            default=lambda o: o.__dict__,
            indent=4,
        )


class RemoteWorkflow:
    """A information container for a remote workflow.

    Args:
        data (dict): workflow information as they are retrieved from the GitHub repository REST API request
            (https://developer.github.com/v3/repos/#get).
    """

    def __init__(self, data):
        # Vars from the initial data payload
        self.name = data.get("name")
        self.full_name = data.get("full_name")
        self.description = data.get("description")
        self.topics = data.get("topics", [])
        self.archived = data.get("archived")
        self.stargazers_count = data.get("stargazers_count")
        self.watchers_count = data.get("watchers_count")
        self.forks_count = data.get("forks_count")

        # Placeholder vars for releases info (ignore pre-releases)
        self.releases = [r for r in data.get("releases", []) if r.get("published_at") is not None]

        # Placeholder vars for local comparison
        self.local_wf = None
        self.local_is_latest = None

        # Beautify date
        for release in self.releases:
            release["published_at_pretty"] = pretty_date(
                datetime.strptime(release.get("published_at"), "%Y-%m-%dT%H:%M:%SZ")
            )
            release["published_at_timestamp"] = int(
                datetime.strptime(release.get("published_at"), "%Y-%m-%dT%H:%M:%SZ").strftime("%s")
            )


class LocalWorkflow:
    """Class to handle local workflows pulled by nextflow"""

    def __init__(self, name):
        """Initialise the LocalWorkflow object"""
        self.full_name = name
        self.repository = None
        self.local_path = None
        self.commit_sha = None
        self.remote_url = None
        self.branch = None
        self.active_tag = None
        self.last_pull = None
        self.last_pull_date = None
        self.last_pull_pretty = None

    def get_local_nf_workflow_details(self):
        """Get full details about a local cached workflow"""

        if self.local_path is None:
            # Try to guess the local cache directory
            if len(os.environ.get("NXF_ASSETS", "")) > 0:
                nf_wfdir = os.path.join(os.environ.get("NXF_ASSETS"), self.full_name)
            elif len(os.environ.get("NXF_HOME", "")) > 0:
                nf_wfdir = os.path.join(os.environ.get("NXF_HOME"), "assets", self.full_name)
            else:
                nf_wfdir = os.path.join(os.getenv("HOME"), ".nextflow", "assets", self.full_name)
            if os.path.isdir(nf_wfdir):
                log.debug(f"Guessed nextflow assets workflow directory: {nf_wfdir}")
                self.local_path = nf_wfdir

            # Use `nextflow info` to get more details about the workflow
            else:
                result = nf_core.utils.run_cmd("nextflow", f"info -d {self.full_name}")
                if result is not None:
                    nfinfo_raw, _ = result
                    re_patterns = {"repository": r"repository\s*: (.*)", "local_path": r"local path\s*: (.*)"}
                    for key, pattern in re_patterns.items():
                        m = re.search(pattern, str(nfinfo_raw))
                        if m:
                            setattr(self, key, m.group(1))

        # Pull information from the local git repository
        if self.local_path is not None:
            log.debug(f"Pulling git info from {self.local_path}")
            try:
                repo = git.Repo(self.local_path)
                self.commit_sha = str(repo.head.commit.hexsha)
                self.remote_url = str(repo.remotes.origin.url)
                self.last_pull = os.stat(os.path.join(self.local_path, ".git", "FETCH_HEAD")).st_mtime
                self.last_pull_date = datetime.fromtimestamp(self.last_pull).strftime("%Y-%m-%d %H:%M:%S")
                self.last_pull_pretty = pretty_date(self.last_pull)

                # Get the checked out branch if we can
                try:
                    self.branch = str(repo.active_branch)
                except TypeError:
                    self.branch = None

                # See if we are on a tag (release)
                self.active_tag = None
                for tag in repo.tags:
                    if str(tag.commit) == str(self.commit_sha):
                        self.active_tag = str(tag)

            # I'm not sure that we need this any more, it predated the self.branch catch above for detacted HEAD
            except (TypeError, git.InvalidGitRepositoryError) as e:
                log.error(
                    f"Could not fetch status of local Nextflow copy of '{self.full_name}':"
                    f"\n   [red]{type(e).__name__}:[/] {str(e)}"
                    "\n\nThis git repository looks broken. It's probably a good idea to delete this local copy and pull again:"
                    f"\n   [magenta]rm -rf {self.local_path}"
                    f"\n   [magenta]nextflow pull {self.full_name}",
                )


def pretty_date(time):
    """Transforms a datetime object or a int() Epoch timestamp into a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc

    Based on https://stackoverflow.com/a/1551394/713980
    Adapted by sven1103
    """

    now = datetime.now()
    if isinstance(time, datetime):
        diff = now - time
    else:
        diff = now - datetime.fromtimestamp(time)
    second_diff = diff.seconds
    day_diff = diff.days

    pretty_msg = (
        (0, ((float("inf"), 1, "from the future"),)),
        (
            1,
            (
                (10, 1, "just now"),
                (60, 1, "{sec:.0f} seconds ago"),
                (120, 1, "a minute ago"),
                (3600, 60, "{sec:.0f} minutes ago"),
                (7200, 1, "an hour ago"),
                (86400, 3600, "{sec:.0f} hours ago"),
            ),
        ),
        (2, ((float("inf"), 1, "yesterday"),)),
        (7, ((float("inf"), 1, "{days:.0f} day{day_s} ago"),)),
        (31, ((float("inf"), 7, "{days:.0f} week{day_s} ago"),)),
        (365, ((float("inf"), 30, "{days:.0f} months ago"),)),
        (float("inf"), ((float("inf"), 365, "{days:.0f} year{day_s} ago"),)),
    )

    for days, seconds in pretty_msg:
        if day_diff < days:
            for sec in seconds:
                if second_diff < sec[0]:
                    return sec[2].format(
                        days=day_diff / sec[1], sec=second_diff / sec[1], day_s="s" if day_diff / sec[1] > 1 else ""
                    )
    return "... time is relative anyway"
