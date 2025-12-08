import logging
import os
import shutil
import tempfile
from pathlib import Path

import questionary

import nf_core.utils
from nf_core.components.components_command import ComponentCommand
from nf_core.components.components_differ import ComponentsDiffer
from nf_core.modules.modules_json import ModulesJson

log = logging.getLogger(__name__)


class ComponentPatch(ComponentCommand):
    def __init__(self, pipeline_dir, component_type, remote_url=None, branch=None, no_pull=False, installed_by=None):
        super().__init__(component_type, pipeline_dir, remote_url, branch, no_pull)

        self.modules_json = ModulesJson(pipeline_dir)

    def _parameter_checks(self, component):
        """Checks the compatibility of the supplied parameters.

        Raises:
            UserWarning: if any checks fail.
        """
        if not self.has_valid_directory():
            raise UserWarning("The command was not run in a valid pipeline directory.")

        components = self.modules_json.get_all_components(self.component_type).get(self.modules_repo.remote_url)
        if components is None:
            raise UserWarning(
                f"No {self.component_type[:-1]}s found in the 'modules.json' file for the remote '{self.modules_repo.remote_url}'"
            )
        component_names = [component for _, component in components]

        if component is not None and component not in component_names:
            component_dir = [d for d, m in components if m == component][0]
            raise UserWarning(
                f"{self.component_type[:-1].title()} '{Path(self.component_type, component_dir, component)}' does not exist in the pipeline"
            )

    def patch(self, component=None):
        # Check modules directory structure
        self.check_modules_structure()

        # Verify that 'modules.json' is consistent with the installed modules
        self.modules_json.check_up_to_date()
        self._parameter_checks(component)
        components = self.modules_json.get_all_components(self.component_type).get(self.modules_repo.remote_url)

        if component is None:
            choices = [
                component if directory == self.modules_repo.repo_path else f"{directory}/{component}"
                for directory, component in components
            ]
            component = questionary.autocomplete(
                f"{self.component_type[:-1].title()} name:",
                choices=sorted(choices),
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
        component_dir = [dir for dir, m in components if m == component][0]
        component_fullname = str(Path(self.component_type, self.modules_repo.repo_path, component))

        # Verify that the component has an entry in the modules.json file
        if not self.modules_json.component_present(
            component, self.modules_repo.remote_url, component_dir, self.component_type
        ):
            raise UserWarning(
                f"The '{component_fullname}' {self.component_type[:-1]} does not have an entry in the 'modules.json' file. Cannot compute patch"
            )

        component_version = self.modules_json.get_component_version(
            self.component_type, component, self.modules_repo.remote_url, self.modules_repo.repo_path
        )
        if component_version is None:
            raise UserWarning(
                f"The '{component_fullname}' {self.component_type[:-1]} does not have a valid version in the 'modules.json' file. Cannot compute patch"
            )
        # Get the component branch and reset it in the ModulesRepo object
        component_branch = self.modules_json.get_component_branch(
            self.component_type, component, self.modules_repo.remote_url, component_dir
        )
        if component_branch != self.modules_repo.branch:
            self.modules_repo.setup_branch(component_branch)

        # Set the diff filename based on the module name
        patch_filename = f"{component.replace('/', '-')}.diff"
        component_relpath = Path(self.component_type, component_dir, component)
        patch_relpath = Path(component_relpath, patch_filename)
        component_current_dir = Path(self.directory, component_relpath)
        patch_path = Path(self.directory, patch_relpath)

        if patch_path.exists():
            remove = questionary.confirm(
                f"Patch exists for {self.component_type[:-1]} '{component_fullname}'. Do you want to regenerate it?",
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
            if remove:
                os.remove(patch_path)
            else:
                return

        # Create a temporary directory for storing the unchanged version of the module
        install_dir = tempfile.mkdtemp()
        component_install_dir = Path(install_dir, component)
        if not self.install_component_files(component, component_version, self.modules_repo, install_dir):
            raise UserWarning(
                f"Failed to install files of {self.component_type[:-1]} '{component}' from remote ({self.modules_repo.remote_url})."
            )

        # Write the patch to a temporary location (otherwise it is printed to the screen later)
        patch_temp_path = tempfile.mktemp()
        try:
            ComponentsDiffer.write_diff_file(
                patch_temp_path,
                component,
                self.modules_repo.repo_path,
                component_install_dir,
                component_current_dir,
                for_git=False,
                dsp_from_dir=component_relpath,
                dsp_to_dir=component_relpath,
            )
            log.debug(f"Patch file wrote to a temporary directory {patch_temp_path}")
        except UserWarning:
            raise UserWarning(f"{self.component_type[:-1]} '{component_fullname}' is unchanged. No patch to compute")

        # Write changes to modules.json
        self.modules_json.add_patch_entry(
            self.component_type, component, self.modules_repo.remote_url, component_dir, patch_relpath
        )
        log.debug(f"Wrote patch path for {self.component_type[:-1]} {component} to modules.json")

        # Show the changes made to the module
        ComponentsDiffer.print_diff(
            component,
            self.modules_repo.repo_path,
            component_install_dir,
            component_current_dir,
            dsp_from_dir=component_current_dir,
            dsp_to_dir=component_current_dir,
        )

        # Finally move the created patch file to its final location
        shutil.move(patch_temp_path, patch_path)
        log.info(f"Patch file of '{component_fullname}' written to '{patch_path}'")

    def remove(self, component):
        # Check modules directory structure
        self.check_modules_structure()

        self.modules_json.check_up_to_date()
        self._parameter_checks(component)
        components = self.modules_json.get_all_components(self.component_type).get(self.modules_repo.remote_url)

        if component is None:
            choices = [
                component if directory == self.modules_repo.repo_path else f"{directory}/{component}"
                for directory, component in components
            ]
            component = questionary.autocomplete(
                f"{self.component_type[:-1]} name:",
                choices,
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
        component_dir = [dir for dir, m in components if m == component][0]
        component_fullname = str(Path(self.component_type, component_dir, component))

        # Verify that the component has an entry in the modules.json file
        if not self.modules_json.component_present(
            component, self.modules_repo.remote_url, component_dir, self.component_type
        ):
            raise UserWarning(
                f"The '{component_fullname}' {self.component_type[:-1]} does not have an entry in the 'modules.json' file. Cannot compute patch"
            )

        component_version = self.modules_json.get_component_version(
            self.component_type, component, self.modules_repo.remote_url, self.modules_repo.repo_path
        )
        if component_version is None:
            raise UserWarning(
                f"The '{component_fullname}' {self.component_type[:-1]} does not have a valid version in the 'modules.json' file. Cannot compute patch"
            )
        # Get the module branch and reset it in the ModulesRepo object
        component_branch = self.modules_json.get_component_branch(
            self.component_type, component, self.modules_repo.remote_url, component_dir
        )
        if component_branch != self.modules_repo.branch:
            self.modules_repo.setup_branch(component_branch)

        # Set the diff filename based on the component name
        patch_filename = f"{component.replace('/', '-')}.diff"
        component_relpath = Path(self.component_type, component_dir, component)
        patch_relpath = Path(component_relpath, patch_filename)
        patch_path = Path(self.directory, patch_relpath)
        component_path = Path(self.directory, component_relpath)

        if patch_path.exists():
            remove = questionary.confirm(
                f"Patch exists for {self.component_type[:-1]} '{component_fullname}'. Are you sure you want to remove?",
                style=nf_core.utils.nfcore_question_style,
            ).unsafe_ask()
            if not remove:
                return

        # Try to apply the patch in reverse and move resulting files to module dir
        temp_component_dir = self.modules_json.try_apply_patch_reverse(
            self.component_type, component, self.modules_repo.repo_path, patch_relpath, component_path
        )
        try:
            for file in Path(temp_component_dir).glob("*"):
                file.rename(component_path.joinpath(file.name))
            os.rmdir(temp_component_dir)
        except Exception as err:
            raise UserWarning(f"There was a problem reverting the patched file: {err}")

        log.info(f"Patch for {component} reverted!")
        # Remove patch file if we could revert the patch
        patch_path.unlink()
        # Write changes to module.json
        self.modules_json.remove_patch_entry(component, self.modules_repo.remote_url, component_dir)

        if not all(
            self.modules_repo.component_files_identical(
                component, component_path, component_version, self.component_type
            ).values()
        ):
            log.error(
                f"Module files do not appear to match the remote for the commit sha in the 'module.json': {component_version}\n"
                f"Recommend reinstalling with 'nf-core modules install --force --sha {component_version} {component}' "
            )
