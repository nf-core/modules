import logging
import os
import shutil
import tempfile
from pathlib import Path

import questionary

import nf_core.modules.modules_utils
import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.components.components_differ import ComponentsDiffer
from nf_core.components.components_utils import (
    get_components_to_install,
    prompt_component_version_sha,
)
from nf_core.components.install import ComponentInstall
from nf_core.components.remove import ComponentRemove
from nf_core.modules.modules_json import ModulesJson
from nf_core.modules.modules_repo import ModulesRepo
from nf_core.utils import plural_es, plural_s, plural_y

log = logging.getLogger(__name__)


class ComponentUpdate(ComponentCommand):
    def __init__(
        self,
        pipeline_dir,
        component_type,
        force=False,
        prompt=False,
        sha=None,
        update_all=False,
        show_diff=None,
        save_diff_fn=None,
        update_deps=False,
        remote_url=None,
        branch=None,
        no_pull=False,
        limit_output=False,
    ):
        super().__init__(component_type, pipeline_dir, remote_url, branch, no_pull)
        self.current_remote = ModulesRepo(remote_url, branch)
        self.branch = branch
        self.force = force
        self.prompt = prompt
        self.sha = sha
        self.update_all = update_all
        self.show_diff = show_diff
        self.save_diff_fn = save_diff_fn
        self.limit_output = limit_output
        self.update_deps = update_deps
        self.component = None
        self.update_config = None
        self.modules_json = ModulesJson(self.directory)
        self.branch = branch

    def _parameter_checks(self):
        """Checks the compatibility of the supplied parameters.

        Raises:
            UserWarning: if any checks fail.
        """

        if self.save_diff_fn and self.show_diff:
            raise UserWarning("Either `--preview` or `--save_diff` can be specified, not both.")

        if self.update_all and self.component:
            raise UserWarning(f"Either a {self.component_type[:-1]} or the '--all' flag can be specified, not both.")

        if self.repo_type == "modules":
            raise UserWarning(
                f"{self.component_type.title()} can not be updated in clones of the nf-core/modules repository."
            )

        if self.prompt and self.sha is not None:
            raise UserWarning("Cannot use '--sha' and '--prompt' at the same time.")

        if not self.has_valid_directory():
            raise UserWarning("The command was not run in a valid pipeline directory.")
        if self.limit_output and not (self.save_diff_fn or self.show_diff):
            raise UserWarning("The '--limit-output' flag can only be used with '--preview' or '--save-diff'.")

    def update(self, component=None, silent=False, updated=None, check_diff_exist=True) -> bool:
        """Updates a specified module/subworkflow or all modules/subworkflows in a pipeline.

        If updating a subworkflow: updates all modules used in that subworkflow.
        If updating a module: updates all subworkflows that use the module.

        Args:
            component (str): The name of the module/subworkflow to update.

        Returns:
            (bool): True if the update was successful, False otherwise.
        """
        if isinstance(component, dict):
            # Override modules_repo when the component to install is a dependency from a subworkflow.
            remote_url = component.get("git_remote", self.current_remote.remote_url)
            branch = component.get("branch", self.branch)
            self.modules_repo = ModulesRepo(remote_url, branch)
            component = component["name"]

        self.component = component
        if updated is None:
            updated = []

        _, tool_config = nf_core.utils.load_tools_config(self.directory)
        self.update_config = getattr(tool_config, "update", {}) or {}
        self._parameter_checks()

        # Check modules directory structure
        self.check_modules_structure()

        # Verify that 'modules.json' is consistent with the installed modules
        if not silent:
            self.modules_json.check_up_to_date()

        if not self.update_all and component is None:
            choices = [f"All {self.component_type}", f"Named {self.component_type[:-1]}"]
            self.update_all = (
                questionary.select(
                    f"Update all {self.component_type} or a single named {self.component_type[:-1]}?",
                    choices=choices,
                    style=nf_core.utils.nfcore_question_style,
                ).unsafe_ask()
                == f"All {self.component_type}"
            )

        # Verify that the provided SHA exists in the repo
        if self.sha is not None and not self.modules_repo.sha_exists_on_branch(self.sha):
            log.error(f"Commit SHA '{self.sha}' doesn't exist in '{self.modules_repo.remote_url}'")
            return False

        # Get the list of modules/subworkflows to update, and their version information
        components_info = (
            self.get_all_components_info() if self.update_all else [self.get_single_component_info(component)]
        )
        # Save the current state of the modules.json
        old_modules_json = self.modules_json.get_modules_json()

        # Ask if we should show the diffs (unless a filename was already given on the command line)
        if not self.save_diff_fn and self.show_diff is None:
            diff_type = questionary.select(
                "Do you want to view diffs of the proposed changes?",
                choices=[
                    {"name": "No previews, just update everything", "value": 0},
                    {"name": "Preview diff in terminal, choose whether to update files", "value": 1},
                    {"name": "Just write diffs to a patch file", "value": 2},
                ],
                default={"name": "No previews, just update everything", "value": 0},
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

            self.show_diff = diff_type == 1
            self.save_diff_fn = diff_type == 2

        if self.save_diff_fn:  # True or a string
            self.setup_diff_file(check_diff_exist)

        # Loop through all components to be updated
        # and do the requested action on them
        exit_value = True
        all_patches_successful = True
        for modules_repo, component, sha, patch_relpath in components_info:
            if component is None:
                # The entry from .nf-core.yml is set to false, skip update of this component
                continue
            component_fullname = str(Path(self.component_type, modules_repo.repo_path, component))
            # Are we updating the files in place or not?
            dry_run = self.show_diff or self.save_diff_fn

            current_version = self.modules_json.get_component_version(
                self.component_type, component, modules_repo.remote_url, modules_repo.repo_path
            )

            # Set the temporary installation folder
            install_tmp_dir = Path(tempfile.mkdtemp())
            component_install_dir = install_tmp_dir / component

            # Compute the component directory
            component_dir = Path(self.directory, self.component_type, modules_repo.repo_path, component)

            if sha is not None:
                version = sha
            elif self.prompt:
                version = prompt_component_version_sha(
                    component, self.component_type, modules_repo=modules_repo, installed_sha=current_version
                )
            else:
                version = modules_repo.get_latest_component_version(component, self.component_type)

            if current_version is not None and not self.force:
                if current_version == version:
                    if self.sha or self.prompt:
                        log.info(f"'{component_fullname}' is already installed at {version}")
                    else:
                        log.info(f"'{component_fullname}' is already up to date")
                    continue

            # Download component files
            if not self.install_component_files(component, version, modules_repo, str(install_tmp_dir)):
                exit_value = False
                continue

            if patch_relpath is not None:
                patch_successful = self.try_apply_patch(
                    component,
                    modules_repo.repo_path,
                    patch_relpath,
                    component_dir,
                    component_install_dir,
                    write_file=False,
                )
                if patch_successful:
                    log.info(f"{self.component_type[:-1].title()} '{component_fullname}' patched successfully")
                else:
                    log.warning(
                        f"Failed to patch {self.component_type[:-1]} '{component_fullname}'. Will proceed with unpatched files."
                    )
                all_patches_successful &= patch_successful

            if dry_run:
                if patch_relpath is not None:
                    if patch_successful:
                        log.info("Current installation is compared against patched version in remote.")
                    else:
                        log.warning("Current installation is compared against unpatched version in remote.")
                # Compute the diffs for the component
                if self.save_diff_fn:
                    log.info(
                        f"Writing diff file for {self.component_type[:-1]} '{component_fullname}' to '{self.save_diff_fn}'"
                    )
                    try:
                        ComponentsDiffer.write_diff_file(
                            self.save_diff_fn,
                            component,
                            modules_repo.repo_path,
                            component_dir,
                            component_install_dir,
                            current_version,
                            version,
                            dsp_from_dir=component_dir,
                            dsp_to_dir=component_dir,
                            limit_output=self.limit_output,
                        )
                        updated.append(component)
                    except UserWarning as e:
                        if str(e) != "Module is unchanged":
                            raise
                        else:
                            updated.append(component)
                    recursive_update = True
                    modules_to_update, subworkflows_to_update = self.get_components_to_update(component)
                    if not silent and len(modules_to_update + subworkflows_to_update) > 0:
                        log.warning(
                            f"All modules and subworkflows linked to the updated {self.component_type[:-1]} will be added to the same diff file.\n"
                            "It is advised to keep all your modules and subworkflows up to date.\n"
                            "It is not guaranteed that a subworkflow will continue working as expected if all modules/subworkflows used in it are not up to date.\n"
                        )
                        if self.update_deps:
                            recursive_update = True
                        else:
                            recursive_update = questionary.confirm(
                                "Would you like to continue adding all modules and subworkflows differences?",
                                default=True,
                                style=nf_core.utils.nfcore_question_style,
                            ).unsafe_ask()
                    if recursive_update and len(modules_to_update + subworkflows_to_update) > 0:
                        # Write all the differences of linked components to a diff file
                        self.update_linked_components(
                            modules_to_update, subworkflows_to_update, updated, check_diff_exist=False
                        )
                        self.manage_changes_in_linked_components(component, modules_to_update, subworkflows_to_update)

                elif self.show_diff:
                    ComponentsDiffer.print_diff(
                        component,
                        modules_repo.repo_path,
                        component_dir,
                        component_install_dir,
                        current_version,
                        version,
                        dsp_from_dir=component_dir,
                        dsp_to_dir=component_dir,
                        limit_output=self.limit_output,
                    )
                    # Ask the user if they want to install the component
                    dry_run = not questionary.confirm(
                        f"Update {self.component_type[:-1]} '{component}'?",
                        default=False,
                        style=nf_core.utils.nfcore_question_style,
                    ).unsafe_ask()

            if not dry_run:
                # Clear the component directory and move the installed files there
                self.move_files_from_tmp_dir(component, str(install_tmp_dir), modules_repo.repo_path, version)
                # Update modules.json with newly installed component
                self.modules_json.update(self.component_type, modules_repo, component, version, installed_by=None)
                updated.append(component)
                recursive_update = True
                modules_to_update, subworkflows_to_update = self.get_components_to_update(component)
                if not silent and len(modules_to_update + subworkflows_to_update) > 0:
                    if not self.update_all:
                        log.warning(
                            f"All modules and subworkflows linked to the updated {self.component_type[:-1]} will be {'asked for update' if self.show_diff else 'automatically updated'}.\n"
                            "It is advised to keep all your modules and subworkflows up to date.\n"
                            "It is not guaranteed that a subworkflow will continue working as expected if all modules/subworkflows used in it are not up to date.\n"
                        )
                        if self.update_deps:
                            recursive_update = True
                        else:
                            recursive_update = questionary.confirm(
                                "Would you like to continue updating all modules and subworkflows?",
                                default=True,
                                style=nf_core.utils.nfcore_question_style,
                            ).unsafe_ask()
                if recursive_update and len(modules_to_update + subworkflows_to_update) > 0:
                    # Update linked components
                    self.update_linked_components(modules_to_update, subworkflows_to_update, updated)
                    self.manage_changes_in_linked_components(component, modules_to_update, subworkflows_to_update)

        if self.save_diff_fn:
            # Write the modules.json diff to the file
            ComponentsDiffer.append_modules_json_diff(
                self.save_diff_fn,
                old_modules_json,
                self.modules_json.get_modules_json(),
                Path(self.directory, "modules.json"),
            )
            if exit_value and not silent:
                log.info(
                    f"[bold magenta italic] TIP! [/] If you are happy with the changes in '{self.save_diff_fn}', you "
                    "can apply them by running the command :point_right:"
                    f"  [bold magenta italic]git apply {self.save_diff_fn} [/]"
                )
        elif not all_patches_successful and not silent:
            log.info(f"Updates complete. Please apply failed patch{plural_es(components_info)} manually.")
            self.modules_json.load()
            self.modules_json.dump(run_prettier=True)
        elif not silent:
            log.info("Updates complete :sparkles:")
            self.modules_json.load()
            self.modules_json.dump(run_prettier=True)

        return exit_value

    def get_single_component_info(self, component):
        """Collects the modules repository, version and sha for a component.

        Information about the component version in the '.nf-core.yml' overrides
        the '--sha' option

        Args:
            component (str): The name of the module/subworkflow to get info for.

        Returns:
            (ModulesRepo, str, str): The modules repo containing the component,
            the component name, and the component version.

        Raises:
            LookupError: If the component is not found either in the pipeline or the modules repo.
            UserWarning: If the '.nf-core.yml' entry is not valid.
        """
        # Check if there are any modules/subworkflows installed from the repo
        repo_url = self.modules_repo.remote_url
        components = self.modules_json.get_all_components(self.component_type).get(repo_url)
        if components is None:
            raise LookupError(f"No {self.component_type} installed from '{repo_url}'")

        choices = [
            component if directory == self.modules_repo.repo_path else f"{directory}/{component}"
            for directory, component in components
        ]

        if component is None:
            component = questionary.autocomplete(
                f"{self.component_type[:-1].title()} name:",
                choices=sorted(choices),
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()

        # Get component installation directory
        try:
            install_dir = [dir for dir, m in components if component == m][0]
        except IndexError:
            raise UserWarning(f"{self.component_type[:-1].title()} '{component}' not found in 'modules.json'.")

        # Check if component is installed before trying to update
        if component not in choices:
            raise LookupError(
                f"{self.component_type[:-1].title()} '{component}' is not installed in pipeline and could therefore not be updated"
            )

        # Check that the supplied name is an available module/subworkflow
        if component and component not in self.modules_repo.get_avail_components(self.component_type, commit=self.sha):
            raise LookupError(
                f"{self.component_type[:-1].title()} '{component}' not found in list of available {self.component_type}."
                f"Use the command 'nf-core {self.component_type} list remote' to view available software"
            )

        sha = self.sha
        config_entry = None
        if self.update_config is not None:
            if any(
                [
                    entry.count("/") == 1
                    and (entry.endswith("modules") or entry.endswith("subworkflows"))
                    and not (entry.endswith(".git") or entry.endswith(".git/"))
                    for entry in self.update_config.keys()
                ]
            ):
                raise UserWarning(
                    "Your '.nf-core.yml' file format is outdated. "
                    "The format should be of the form:\n"
                    "update:\n  <repo_url>:\n    <component_install_directory>:\n      <component_name>:"
                )
            if isinstance(self.update_config.get(self.modules_repo.remote_url, {}), str):
                # If the repo entry is a string, it's the sha to update to
                config_entry = self.update_config.get(self.modules_repo.remote_url, {})
            elif component in self.update_config.get(self.modules_repo.remote_url, {}).get(install_dir, {}):
                # If the component to update is in .nf-core.yml config file
                config_entry = self.update_config[self.modules_repo.remote_url][install_dir].get(component)
        if config_entry is not None and config_entry is not True:
            if config_entry is False:
                log.warn(
                    f"{self.component_type[:-1].title()}'s update entry in '.nf-core.yml' for '{component}' is set to False"
                )
                return (self.modules_repo, None, None, None)
            if not isinstance(config_entry, str):
                raise UserWarning(
                    f"{self.component_type[:-1].title()}'s update entry in '.nf-core.yml' for '{component}' is of wrong type"
                )

            sha = config_entry
            if self.sha is not None:
                log.warning(
                    f"Found entry in '.nf-core.yml' for {self.component_type[:-1]} '{component}' "
                    "which will override version specified with '--sha'"
                )
            else:
                log.info(f"Found entry in '.nf-core.yml' for {self.component_type[:-1]} '{component}'")
            log.info(f"Updating component to ({sha})")

        # Check if the update branch is the same as the installation branch
        current_branch = self.modules_json.get_component_branch(
            self.component_type, component, self.modules_repo.remote_url, install_dir
        )
        new_branch = self.modules_repo.branch
        if current_branch != new_branch:
            log.warning(
                f"You are trying to update the '{Path(install_dir, component)}' {self.component_type[:-1]} from "
                f"the '{new_branch}' branch. This {self.component_type[:-1]} was installed from the '{current_branch}'"
            )
            switch = questionary.confirm(f"Do you want to update using the '{current_branch}' instead?").unsafe_ask()
            if switch:
                # Change the branch
                self.modules_repo.setup_branch(current_branch)

        # If there is a patch file, get its filename
        patch_fn = self.modules_json.get_patch_fn(
            self.component_type, component, self.modules_repo.remote_url, install_dir
        )

        return (self.modules_repo, component, sha, patch_fn)

    def get_all_components_info(self, branch=None):
        """Collects the modules repository, version and sha for all modules/subworkflows.

        Information about the module/subworkflow version in the '.nf-core.yml' overrides the '--sha' option.

        Returns:
            [(ModulesRepo, str, str)]: A list of tuples containing a ModulesRepo object,
            the component name, and the component version.
        """
        if branch is not None:
            use_branch = questionary.confirm(
                f"'--branch' was specified. Should this branch be used to update all {self.component_type}?",
                default=False,
            )
            if not use_branch:
                branch = None
        skipped_repos = []
        skipped_components = []
        overridden_repos = []
        overridden_components = []
        components_info = {}
        # Loop through all the modules/subworkflows in the pipeline
        # and check if they have an entry in the '.nf-core.yml' file

        for repo_name, components in self.modules_json.get_all_components(self.component_type).items():
            if isinstance(self.update_config, dict) and (
                repo_name not in self.update_config or self.update_config[repo_name] is True
            ):
                # There aren't restrictions for the repository in .nf-core.yml file
                components_info[repo_name] = {}
                for component_dir, component in components:
                    try:
                        components_info[repo_name][component_dir].append(
                            (
                                component,
                                self.sha,
                                self.modules_json.get_component_branch(
                                    self.component_type, component, repo_name, component_dir
                                ),
                            )
                        )
                    except KeyError:
                        components_info[repo_name][component_dir] = [
                            (
                                component,
                                self.sha,
                                self.modules_json.get_component_branch(
                                    self.component_type, component, repo_name, component_dir
                                ),
                            )
                        ]
            elif isinstance(self.update_config, dict) and isinstance(self.update_config[repo_name], dict):
                # If it is a dict, then there are entries for individual components or component directories
                for component_dir in set([dir for dir, _ in components]):
                    if isinstance(self.update_config[repo_name][component_dir], str):
                        # If a string is given it is the commit SHA to which we should update to
                        custom_sha = self.update_config[repo_name][component_dir]
                        components_info[repo_name] = {}
                        for dir, component in components:
                            if component_dir == dir:
                                try:
                                    components_info[repo_name][component_dir].append(
                                        (
                                            component,
                                            custom_sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    )
                                except KeyError:
                                    components_info[repo_name][component_dir] = [
                                        (
                                            component,
                                            custom_sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    ]
                        if self.sha is not None:
                            overridden_repos.append(repo_name)
                    elif self.update_config[repo_name][component_dir] is False:
                        for directory, component in components:
                            if directory == component_dir:
                                skipped_components.append(f"{component_dir}/{components}")
                    elif isinstance(self.update_config[repo_name][component_dir], dict):
                        # If it's a dict, there are entries for individual components
                        dir_config = self.update_config[repo_name][component_dir]
                        components_info[repo_name] = {}
                        for component_dir, component in components:
                            if component not in dir_config or dir_config[component] is True:
                                try:
                                    components_info[repo_name][component_dir].append(
                                        (
                                            component,
                                            self.sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    )
                                except KeyError:
                                    components_info[repo_name][component_dir] = [
                                        (
                                            component,
                                            self.sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    ]
                            elif isinstance(dir_config[component], str):
                                # If a string is given it is the commit SHA to which we should update to
                                custom_sha = dir_config[component]
                                try:
                                    components_info[repo_name][component_dir].append(
                                        (
                                            component,
                                            custom_sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    )
                                except KeyError:
                                    components_info[repo_name][component_dir] = [
                                        (
                                            component,
                                            custom_sha,
                                            self.modules_json.get_component_branch(
                                                self.component_type, component, repo_name, component_dir
                                            ),
                                        )
                                    ]
                                if self.sha is not None:
                                    overridden_components.append(component)
                            elif dir_config[component] is False:
                                # Otherwise the entry must be 'False' and we should ignore the component
                                skipped_components.append(f"{component_dir}/{component}")
                            else:
                                raise UserWarning(
                                    f"{self.component_type[:-1].title()} '{component}' in '{component_dir}' has an invalid entry in '.nf-core.yml'"
                                )
            elif isinstance(self.update_config, dict) and isinstance(self.update_config[repo_name], str):
                # If a string is given it is the commit SHA to which we should update to
                custom_sha = self.update_config[repo_name]
                components_info[repo_name] = {}
                for component_dir, component in components:
                    try:
                        components_info[repo_name][component_dir].append(
                            (
                                component,
                                custom_sha,
                                self.modules_json.get_component_branch(
                                    self.component_type, component, repo_name, component_dir
                                ),
                            )
                        )
                    except KeyError:
                        components_info[repo_name][component_dir] = [
                            (
                                component,
                                custom_sha,
                                self.modules_json.get_component_branch(
                                    self.component_type, component, repo_name, component_dir
                                ),
                            )
                        ]
                if self.sha is not None:
                    overridden_repos.append(repo_name)
            elif isinstance(self.update_config, dict) and self.update_config[repo_name] is False:
                skipped_repos.append(repo_name)
            else:
                log.debug(f"no update config for {repo_name} in `.nf-core.yml`")

        if skipped_repos:
            skipped_str = "', '".join(skipped_repos)
            log.info(f"Skipping {self.component_type} in repositor{plural_y(skipped_repos)}: '{skipped_str}'")

        if skipped_components:
            skipped_str = "', '".join(skipped_components)
            log.info(f"Skipping {self.component_type[:-1]}{plural_s(skipped_components)}: '{skipped_str}'")

        if overridden_repos:
            overridden_str = "', '".join(overridden_repos)
            log.info(
                f"Overriding '--sha' flag for {self.component_type} in repositor{plural_y(overridden_repos)} "
                f"with '.nf-core.yml' entry: '{overridden_str}'"
            )

        if overridden_components:
            overridden_str = "', '".join(overridden_components)
            log.info(
                f"Overriding '--sha' flag for {self.component_type[:-1]}{plural_s(overridden_components)} with "
                f"'.nf-core.yml' entry: '{overridden_str}'"
            )
        # Loop through components_info and create on ModulesRepo object per remote and branch
        repos_and_branches = {}
        for repo_name, repo_content in components_info.items():
            for component_dir, comps in repo_content.items():
                for comp, sha, comp_branch in comps:
                    if branch is not None:
                        comp_branch = branch
                    if (repo_name, comp_branch) not in repos_and_branches:
                        repos_and_branches[(repo_name, comp_branch)] = []
                    repos_and_branches[(repo_name, comp_branch)].append((comp, sha))

        # Create ModulesRepo objects
        repo_objs_comps = []
        for (repo_url, branch), comps_shas in repos_and_branches.items():
            try:
                modules_repo = ModulesRepo(remote_url=repo_url, branch=branch)
            except LookupError as e:
                log.warning(e)
                log.info(f"Skipping {self.component_type} in '{repo_url}'")
            else:
                repo_objs_comps.append((modules_repo, comps_shas))

        # Flatten the list
        components_info = [(repo, comp, sha) for repo, comps_shas in repo_objs_comps for comp, sha in comps_shas]

        # Verify that that all components and shas exist in their respective ModulesRepo,
        # don't try to update those that don't
        i = 0
        while i < len(components_info):
            repo, component, sha = components_info[i]
            if not repo.component_exists(component, self.component_type):
                log.warning(
                    f"{self.component_type[:-1].title()} '{component}' does not exist in '{repo.remote_url}'. Skipping..."
                )
                components_info.pop(i)
            elif sha is not None and not repo.sha_exists_on_branch(sha):
                log.warning(
                    f"Git sha '{sha}' does not exists on the '{repo.branch}' of '{repo.remote_url}'. Skipping {self.component_type[:-1]} '{component}'"
                )
                components_info.pop(i)
            else:
                i += 1

        # Add patch filenames to the components that have them
        components_info = [
            (
                repo,
                comp,
                sha,
                self.modules_json.get_patch_fn(self.component_type, comp, repo.remote_url, repo.repo_path),
            )
            for repo, comp, sha in components_info
        ]

        return components_info

    def setup_diff_file(self, check_diff_exist=True):
        """Sets up the diff file.

        If the save diff option was chosen interactively, the user is asked to supply a name for the diff file.

        Then creates the file for saving the diff.
        """
        if self.save_diff_fn is True:
            # From questionary - no filename yet
            self.save_diff_fn = questionary.path(
                "Enter the filename: ", style=nf_core.utils.nfcore_question_style
            ).unsafe_ask()
        if self.save_diff_fn is not None:
            self.save_diff_fn = Path(self.save_diff_fn)
        else:
            raise UserWarning("No filename provided for saving the diff file")

        if not check_diff_exist:
            # This guarantees that the file exists after calling the function
            self.save_diff_fn.touch()
            return
        # Check if filename already exists (questionary or cli)
        while self.save_diff_fn.exists():
            if questionary.confirm(f"'{self.save_diff_fn}' exists. Remove file?").unsafe_ask():
                os.remove(self.save_diff_fn)
                break
            self.save_diff_fn = questionary.path(
                "Enter a new filename: ",
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
            self.save_diff_fn = Path(self.save_diff_fn)

        # This guarantees that the file exists after calling the function
        self.save_diff_fn.touch()

    def move_files_from_tmp_dir(self, component: str, install_folder: str, repo_path: str, new_version: str) -> None:
        """Move the files from the temporary to the installation directory.

        Args:
            component (str): The module/subworkflow name.
            install_folder [str]: The path to the temporary installation directory.
            repo_path (str): The name of the directory where modules/subworkflows are installed
            new_version (str): The version of the module/subworkflow that was installed.
        """
        temp_component_dir = Path(install_folder, component)
        files = [file_path for file_path in temp_component_dir.rglob("*") if file_path.is_file()]
        pipeline_path = Path(self.directory, self.component_type, repo_path, component)

        if pipeline_path.exists():
            pipeline_files = [f.name for f in pipeline_path.iterdir() if f.is_file()]
            # check if any *.config file exists in the pipeline
            config_files = [f for f in pipeline_files if str(f).endswith(".config")]
            for config_file in config_files:
                log.debug(f"Moving '{component}/{config_file}' to updated component")
                shutil.move(str(pipeline_path / config_file), temp_component_dir / config_file)
                files.append(temp_component_dir / config_file)

        else:
            log.debug(f"Creating new {self.component_type[:-1]} '{component}' in '{self.component_type}/{repo_path}'")

        log.debug(f"Removing old version of {self.component_type[:-1]} '{component}'")
        self.clear_component_dir(component, str(pipeline_path))

        pipeline_path.mkdir(parents=True, exist_ok=True)
        for file in files:
            file = file.relative_to(temp_component_dir)
            path = Path(temp_component_dir, file)
            if path.exists():
                log.debug(f"Moving '{file}' to updated component")
                dest = Path(pipeline_path, file)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(path), dest)

        log.info(f"Updating '{repo_path}/{component}'")
        log.debug(f"Updating {self.component_type[:-1]} '{component}' to {new_version} from {repo_path}")

    def try_apply_patch(
        self, component, repo_path, patch_relpath, component_dir, component_install_dir, write_file=True
    ):
        """
        Try applying a patch file to the new module/subworkflow files


        Args:
            component (str): The name of the module/subworkflow
            repo_path (str): The name of the repository where the module/subworkflow resides
            patch_relpath (Path | str): The path to patch file in the pipeline
            component_dir (Path | str): The module/subworkflow directory in the pipeline
            component_install_dir (Path | str): The directory where the new component
                                            file have been installed

        Returns:
            (bool): Whether the patch application was successful
        """
        component_fullname = str(Path(repo_path, component))
        log.info(f"Found patch for  {self.component_type[:-1]} '{component_fullname}'. Trying to apply it to new files")

        patch_path = Path(self.directory / patch_relpath)
        component_relpath = Path(self.component_type, repo_path, component)

        # Check that paths in patch file are updated
        self.check_patch_paths(patch_path, component)

        # Copy the installed files to a new temporary directory to save them for later use
        temp_dir = Path(tempfile.mkdtemp())
        temp_component_dir = temp_dir / component
        shutil.copytree(component_install_dir, temp_component_dir)

        try:
            new_files = ComponentsDiffer.try_apply_patch(
                self.component_type, component, repo_path, patch_path, temp_component_dir
            )
        except LookupError:
            # Patch failed. Save the patch file by moving to the install dir
            shutil.move(patch_path, Path(component_install_dir, patch_path.relative_to(component_dir)))
            log.warning(
                f"Failed to apply patch for {self.component_type[:-1]} '{component_fullname}'. You will have to apply the patch manually"
            )
            return False

        # Write the patched files to a temporary directory
        log.debug("Writing patched files")
        for file, new_content in new_files.items():
            fn = temp_component_dir / file
            with open(fn, "w") as fh:
                fh.writelines(new_content)

        # Create the new patch file
        log.debug("Regenerating patch file")
        ComponentsDiffer.write_diff_file(
            Path(temp_component_dir, patch_path.relative_to(component_dir)),
            component,
            repo_path,
            component_install_dir,
            temp_component_dir,
            file_action="w",
            for_git=False,
            dsp_from_dir=component_relpath,
            dsp_to_dir=component_relpath,
            limit_output=self.limit_output,
        )

        # Move the patched files to the install dir
        log.debug("Overwriting installed files installed files  with patched files")
        shutil.rmtree(component_install_dir)
        shutil.copytree(temp_component_dir, component_install_dir)

        # Add the patch file to the modules.json file
        self.modules_json.add_patch_entry(
            self.component_type,
            component,
            self.modules_repo.remote_url,
            repo_path,
            patch_relpath,
            write_file=write_file,
        )

        return True

    def get_components_to_update(self, component):
        """
        Get all modules and subworkflows linked to the updated component.

        Returns:
            (list,list): A tuple of lists with the modules and subworkflows to update
        """
        mods_json = self.modules_json.get_modules_json()
        modules_to_update = []
        subworkflows_to_update = []
        installed_by = self.modules_json.get_installed_by_entries(self.component_type, component)

        if self.component_type == "modules":
            # All subworkflow names in the installed_by section of a module are subworkflows using this module
            # We need to update them too
            git_remote = self.current_remote.remote_url
            for subworkflow in installed_by:
                if subworkflow != component:
                    for remote_url, content in mods_json["repos"].items():
                        if all_subworkflows := content.get("subworkflows"):
                            for details in all_subworkflows.values():
                                if subworkflow in details:
                                    git_remote = remote_url
                if subworkflow != self.component_type:
                    subworkflows_to_update.append({"name": subworkflow, "git_remote": git_remote})

        elif self.component_type == "subworkflows":
            for repo, repo_content in mods_json["repos"].items():
                for component_type, dir_content in repo_content.items():
                    for dir, components in dir_content.items():
                        for comp, comp_content in components.items():
                            # If the updated subworkflow name appears in the installed_by section of the checked component
                            # The checked component is used by the updated subworkflow
                            # We need to update it too
                            if component in comp_content["installed_by"]:
                                if component_type == "modules":
                                    modules_to_update.append({"name": comp, "git_remote": repo, "org_path": dir})
                                elif component_type == "subworkflows":
                                    subworkflows_to_update.append({"name": comp, "git_remote": repo, "org_path": dir})

        return modules_to_update, subworkflows_to_update

    def update_linked_components(
        self,
        modules_to_update,
        subworkflows_to_update,
        updated=None,
        check_diff_exist=True,
    ):
        """
        Update modules and subworkflows linked to the component being updated.
        """
        for s_update in subworkflows_to_update:
            if s_update["name"] in updated:
                continue
            original_component_type, original_update_all = self._change_component_type("subworkflows")
            self.update(
                s_update,
                silent=True,
                updated=updated,
                check_diff_exist=check_diff_exist,
            )
            self._reset_component_type(original_component_type, original_update_all)

        for m_update in modules_to_update:
            if m_update["name"] in updated:
                continue
            original_component_type, original_update_all = self._change_component_type("modules")
            try:
                self.update(
                    m_update,
                    silent=True,
                    updated=updated,
                    check_diff_exist=check_diff_exist,
                )
            except LookupError as e:
                # If the module to be updated is not available, check if there has been a name change
                if "not found in list of available" in str(e):
                    # Skip update, we check for name changes with manage_changes_in_linked_components
                    pass
                else:
                    raise
            finally:
                self._reset_component_type(original_component_type, original_update_all)

    def manage_changes_in_linked_components(self, component, modules_to_update, subworkflows_to_update):
        """Check for linked components added or removed in the new subworkflow version"""
        if self.component_type == "subworkflows":
            org_path = self.current_remote.repo_path

            subworkflow_directory = Path(self.directory, self.component_type, org_path, component)
            included_modules, included_subworkflows = get_components_to_install(subworkflow_directory)
            # If a module/subworkflow has been removed from the subworkflow
            for module in modules_to_update:
                module_name = module["name"]
                included_modules_names = [m["name"] for m in included_modules]
                if module_name not in included_modules_names:
                    log.info(f"Removing module '{module_name}' which is not included in '{component}' anymore.")
                    remove_module_object = ComponentRemove("modules", self.directory)
                    remove_module_object.remove(module_name, removed_by=component)
            for subworkflow in subworkflows_to_update:
                subworkflow_name = subworkflow["name"]
                included_subworkflow_names = [m["name"] for m in included_subworkflows]
                if subworkflow_name not in included_subworkflow_names:
                    log.info(
                        f"Removing subworkflow '{subworkflow_name}' which is not included in '{component}' anymore."
                    )
                    remove_subworkflow_object = ComponentRemove("subworkflows", self.directory)
                    remove_subworkflow_object.remove(subworkflow_name, removed_by=component)
            # If a new module/subworkflow is included in the subworklfow and wasn't included before
            for module in included_modules:
                module_name = module["name"]
                module["git_remote"] = module.get("git_remote", self.current_remote.remote_url)
                module["branch"] = module.get("branch", self.branch)
                if module_name not in modules_to_update:
                    log.info(f"Installing newly included module '{module_name}' for '{component}'")
                    install_module_object = ComponentInstall(self.directory, "modules", installed_by=component)
                    install_module_object.install(module, silent=True)
            for subworkflow in included_subworkflows:
                subworkflow_name = subworkflow["name"]
                subworkflow["git_remote"] = subworkflow.get("git_remote", self.current_remote.remote_url)
                subworkflow["branch"] = subworkflow.get("branch", self.branch)
                if subworkflow_name not in subworkflows_to_update:
                    log.info(f"Installing newly included subworkflow '{subworkflow_name}' for '{component}'")
                    install_subworkflow_object = ComponentInstall(
                        self.directory, "subworkflows", installed_by=component
                    )
                    install_subworkflow_object.install(subworkflow, silent=True)

    def _change_component_type(self, new_component_type):
        original_component_type = self.component_type
        self.component_type = new_component_type
        self.modules_json.pipeline_components = None
        # also reset update_all in case it's set
        original_update_all = self.update_all
        self.update_all = False
        return original_component_type, original_update_all

    def _reset_component_type(self, original_component_type, original_update_all):
        self.component_type = original_component_type
        self.modules_json.pipeline_components = None
        self.update_all = original_update_all
        if self.current_remote is None:
            self.current_remote = self.modules_repo
