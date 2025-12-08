import json
import logging
from pathlib import Path
from typing import cast

import rich.table

from nf_core.components.components_command import ComponentCommand
from nf_core.modules.modules_json import ModulesJson, ModulesJsonModuleEntry
from nf_core.modules.modules_repo import ModulesRepo

log = logging.getLogger(__name__)


class ComponentList(ComponentCommand):
    def __init__(
        self,
        component_type: str,
        pipeline_dir: str | Path = ".",
        remote: bool = True,
        remote_url: str | None = None,
        branch: str | None = None,
        no_pull: bool = False,
    ) -> None:
        self.remote = remote
        super().__init__(component_type, pipeline_dir, remote_url, branch, no_pull)

    def _configure_repo_and_paths(self, nf_dir_req: bool = True) -> None:
        """
        Override the default with nf_dir_req set to False to allow
        info to be run from anywhere and still return remote info
        """
        if self.remote:
            nf_dir_req = False
        return super()._configure_repo_and_paths(nf_dir_req)

    def list_components(self, keywords: list[str] | None = None, print_json: bool = False) -> rich.table.Table | str:
        keywords = keywords or []
        """
        Get available modules/subworkflows names from GitHub tree for repo
        and print as list to stdout
        """
        # Check modules directory structure
        # self.check_component_structure(self.component_type)

        # Initialise rich table
        table: rich.table.Table = rich.table.Table()
        table.add_column(f"{self.component_type[:-1].capitalize()} Name")
        components: list[str] = []

        def pattern_msg(keywords: list[str]) -> str:
            if len(keywords) == 0:
                return ""
            if len(keywords) == 1:
                return f" matching pattern '{keywords[0]}'"
            else:
                quoted_keywords = (f"'{key}'" for key in keywords)
                return f" matching patterns {', '.join(quoted_keywords)}"

        # No pipeline given - show all remote
        if self.remote:
            # Filter the modules/subworkflows by keywords
            components = [
                comp
                for comp in self.modules_repo.get_avail_components(self.component_type)
                if all(k in comp for k in keywords)
            ]

            # Nothing found
            if len(components) == 0:
                log.info(
                    f"No available {self.component_type} found in {self.modules_repo.remote_url} ({self.modules_repo.branch})"
                    f"{pattern_msg(keywords)}"
                )
                return ""

            for comp in sorted(components):
                table.add_row(comp)

        # We have a pipeline - list what's installed
        else:
            # Check that we are in a pipeline directory
            log.info(f"Repository type: [blue]{self.repo_type}")
            try:
                if self.repo_type != "pipeline":
                    raise UserWarning(
                        f"The command 'nf-core {self.component_type} list local' must be run from a pipeline directory.",
                    )
            except UserWarning as e:
                log.error(e)
                return ""
            # Check whether pipelines is valid
            try:
                self.has_valid_directory()
            except UserWarning as e:
                log.error(e)
                return ""

            # Verify that 'modules.json' is consistent with the installed modules
            modules_json: ModulesJson = ModulesJson(self.directory)
            modules_json.check_up_to_date()

            # Filter by keywords
            repos_with_comps = {
                repo_url: [comp for comp in components if all(k in comp[1] for k in keywords)]
                for repo_url, components in modules_json.get_all_components(self.component_type).items()
            }

            # Nothing found
            if sum(map(len, repos_with_comps)) == 0:
                log.info(f"No nf-core {self.component_type} found in '{self.directory}'{pattern_msg(keywords)}")
                return ""

            table.add_column("Repository")
            table.add_column("Version SHA")
            table.add_column("Message")
            table.add_column("Date")

            # Load 'modules.json'
            modules_json_file = modules_json.modules_json

            for repo_url, component_with_dir in sorted(repos_with_comps.items()):
                repo_entry: dict[str, dict[str, dict[str, ModulesJsonModuleEntry]]]
                if modules_json_file is None:
                    log.warning(f"Modules JSON file '{modules_json.modules_json_path}' is missing. ")
                    continue
                else:
                    repo_entry = modules_json_file["repos"].get(repo_url, {})
                    for install_dir, component in sorted(component_with_dir):
                        # Use cast() to predict the return type of recursive get():s
                        repo_modules = cast(dict, repo_entry.get(self.component_type))
                        component_entry = cast(dict, cast(dict, repo_modules.get(install_dir)).get(component))

                        if component_entry:
                            version_sha = component_entry["git_sha"]
                            try:
                                # pass repo_name to get info on modules even outside nf-core/modules
                                module = ModulesRepo(
                                    remote_url=repo_url,
                                    branch=component_entry["branch"],
                                )
                                message, date = module.get_commit_info(version_sha)
                            except LookupError as e:
                                log.warning(e)
                                date = "[red]Not Available"
                                message = "[red]Not Available"
                        else:
                            log.warning(
                                f"Commit SHA for {self.component_type[:-1]} '{install_dir}/{self.component_type}' is missing from 'modules.json'"
                            )
                            version_sha = "[red]Not Available"
                            date = "[red]Not Available"
                            message = "[red]Not Available"
                        nice_repo_name = repo_url.replace("https://github.com/", "").replace(".git", "")
                        table.add_row(
                            component,
                            f"[link={repo_url}]{nice_repo_name}[/link]",
                            f"[link={module.gitless_repo()}/commit/{version_sha}]{version_sha[:7]}[/link]",
                            message,
                            date,
                        )
                        components.append(component)

        if print_json:
            return json.dumps(components, sort_keys=True, indent=4)

        if self.remote:
            log.info(
                f"{self.component_type.capitalize()} available from {self.modules_repo.remote_url} ({self.modules_repo.branch})"
                f"{pattern_msg(keywords)}:\n"
            )
        else:
            log.info(f"{self.component_type.capitalize()} installed in '{self.directory}'{pattern_msg(keywords)}:\n")
        return table
