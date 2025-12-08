import io

import pytest
import ruamel.yaml

import nf_core.modules.lint
from nf_core.components.lint import ComponentLint
from nf_core.components.nfcore_component import NFCoreComponent
from nf_core.modules.lint.environment_yml import environment_yml

from ...test_modules import TestModules

yaml = ruamel.yaml.YAML()
yaml.indent(mapping=2, sequence=2, offset=2)


def yaml_dump_to_string(data):
    stream = io.StringIO()
    yaml.dump(data, stream)
    return stream.getvalue()


class DummyModule(NFCoreComponent):
    def __init__(self, path):
        self.environment_yml = path
        self.component_dir = path.parent
        self.component_name = "dummy"
        self.passed = []
        self.failed = []
        self.warned = []


class DummyLint(ComponentLint):
    def __init__(self, tmp_path):
        self.modules_repo = type("repo", (), {"local_repo_dir": tmp_path})
        self.passed = []
        self.failed = []


def setup_test_environment(tmp_path, content, filename="environment.yml"):
    test_file = tmp_path / filename
    test_file.write_text(content)

    (tmp_path / "modules").mkdir(exist_ok=True)
    (tmp_path / "modules" / "environment-schema.json").write_text("{}")

    module = DummyModule(test_file)
    lint = DummyLint(tmp_path)

    return test_file, module, lint


def assert_yaml_result(test_file, expected):
    result = test_file.read_text()
    lines = result.splitlines(True)

    if lines[:2] == [
        "---\n",
        "# yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json\n",
    ]:
        parsed = yaml.load("".join(lines[2:]))
    else:
        parsed = yaml.load(result)

    if isinstance(expected, list):
        assert parsed["dependencies"] == expected
    else:
        for key, value in expected.items():
            assert key in parsed
            assert parsed[key] == value


@pytest.mark.parametrize(
    "input_content,expected",
    [
        # Test basic dependency sorting
        (
            """
            dependencies:
                - zlib
                - python
            """,
            ["python", "zlib"],
        ),
        # Test dict dependency sorting
        (
            """
            dependencies:
                - pip:
                    - b
                    - a
                - python
            """,
            ["python", {"pip": ["a", "b"]}],
        ),
        # Test existing headers
        ("---\n# yaml-language-server: $schema=...\ndependencies:\n  - b\n  - a\n", ["a", "b"]),
        # Test channel preservation (no sorting) - channels order preserved as per nf-core/modules#8554
        (
            """
            channels:
                - conda-forge
                - bioconda
            dependencies:
                - python
            """,
            {
                "channels": ["conda-forge", "bioconda"],
                "dependencies": ["python"],
            },
        ),
        # Test channel preservation with additional channels - channels order preserved as per nf-core/modules#8554
        (
            """
            channels:
                - bioconda
                - conda-forge
                - defaults
                - r
            """,
            {
                "channels": ["bioconda", "conda-forge", "defaults", "r"],
            },
        ),
        # Test namespaced dependencies
        (
            """
            dependencies:
                - bioconda::ngscheckmate=1.0.1
                - bioconda::bcftools=1.21
            """,
            ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1"],
        ),
        # Test mixed dependencies
        (
            """
            dependencies:
                - bioconda::ngscheckmate=1.0.1
                - python
                - bioconda::bcftools=1.21
            """,
            ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1", "python"],
        ),
        # Test full environment with channels and namespaced dependencies - channels order preserved as per nf-core/modules#8554
        (
            """
            channels:
                - conda-forge
                - bioconda
            dependencies:
                - bioconda::ngscheckmate=1.0.1
                - bioconda::bcftools=1.21
            """,
            {
                "channels": ["conda-forge", "bioconda"],
                "dependencies": ["bioconda::bcftools=1.21", "bioconda::ngscheckmate=1.0.1"],
            },
        ),
    ],
)
def test_environment_yml_sorting(tmp_path, input_content, expected):
    """Test that environment.yml files are sorted correctly"""
    test_file, module, lint = setup_test_environment(tmp_path, input_content)

    environment_yml(lint, module)

    assert_yaml_result(test_file, expected)
    assert any("environment_yml_sorted" in x for x in [p[1] for p in module.passed])


@pytest.mark.parametrize(
    "invalid_content,filename",
    [
        ("invalid: yaml: here", "bad.yml"),
        ("", "empty.yml"),
    ],
)
def test_environment_yml_invalid_files(tmp_path, invalid_content, filename):
    """Test that invalid YAML files raise exceptions"""
    test_file, module, lint = setup_test_environment(tmp_path, invalid_content, filename)

    with pytest.raises(Exception):
        environment_yml(lint, module)


def test_environment_yml_missing_dependencies(tmp_path):
    """Test handling of environment.yml without dependencies section"""
    content = "channels:\n  - conda-forge\n"
    test_file, module, lint = setup_test_environment(tmp_path, content)

    environment_yml(lint, module)

    expected = {"channels": ["conda-forge"]}
    assert_yaml_result(test_file, expected)


# Integration tests using the full ModuleLint class
@pytest.mark.integration
class TestModulesEnvironmentYmlIntegration(TestModules):
    """Integration tests for environment.yml linting using real modules"""

    def test_modules_environment_yml_file_doesnt_exists(self):
        """Test linting a module with an environment.yml file"""
        # Use context manager for file manipulation
        backup_path = self.bpipe_test_module_path / "environment.yml.bak"
        env_path = self.bpipe_test_module_path / "environment.yml"

        env_path.rename(backup_path)
        try:
            module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
            module_lint.lint(print_results=False, module="bpipe/test")

            assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
            assert len(module_lint.passed) > 0
            assert len(module_lint.warned) >= 0
            assert module_lint.failed[0].lint_test == "environment_yml_exists"
        finally:
            backup_path.rename(env_path)

    def test_modules_environment_yml_file_sorted_correctly(self):
        """Test linting a module with a correctly sorted environment.yml file"""
        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_environment_yml_file_sorted_incorrectly(self):
        """Test linting a module with an incorrectly sorted environment.yml file"""
        env_path = self.bpipe_test_module_path / "environment.yml"

        # Read, modify, and write back
        with open(env_path) as fh:
            yaml_content = yaml.load(fh)

        # Add a new dependency to the environment.yml file and reverse the order
        yaml_content["dependencies"].append("z=0.0.0")
        yaml_content["dependencies"].reverse()

        with open(env_path, "w") as fh:
            fh.write(yaml_dump_to_string(yaml_content))

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        # we fix the sorting on the fly, so this should pass
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_environment_yml_file_dependencies_not_array(self):
        """Test linting a module with dependencies not as an array"""
        env_path = self.bpipe_test_module_path / "environment.yml"

        with open(env_path) as fh:
            yaml_content = yaml.load(fh)

        yaml_content["dependencies"] = "z"

        with open(env_path, "w") as fh:
            fh.write(yaml_dump_to_string(yaml_content))

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "environment_yml_valid"

    def test_modules_environment_yml_file_mixed_dependencies(self):
        """Test linting a module with mixed-type dependencies (strings and pip dict)"""
        env_path = self.bpipe_test_module_path / "environment.yml"

        with open(env_path) as fh:
            yaml_content = yaml.load(fh)

        # Create mixed dependencies with strings and pip dict in wrong order
        yaml_content["dependencies"] = [
            "python=3.8",
            "bioconda::samtools=1.15.1",
            "bioconda::fastqc=0.12.1",
            "pip=23.3.1",
            {"pip": ["zzz-package==1.0.0", "aaa-package==2.0.0"]},
        ]

        with open(env_path, "w") as fh:
            fh.write(yaml_dump_to_string(yaml_content))

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        # Check that the dependencies were sorted correctly
        with open(env_path) as fh:
            sorted_yaml = yaml.load(fh)

        expected_deps = [
            "bioconda::fastqc=0.12.1",
            "bioconda::samtools=1.15.1",
            "python=3.8",
            "pip=23.3.1",
            {"pip": ["aaa-package==2.0.0", "zzz-package==1.0.0"]},
        ]

        assert sorted_yaml["dependencies"] == expected_deps
        assert len(module_lint.failed) == 0, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0

    def test_modules_environment_yml_file_default_channel_fails(self):
        """Test linting a module with invalid default channel in the environment.yml file"""
        env_path = self.bpipe_test_module_path / "environment.yml"

        with open(env_path) as fh:
            yaml_content = yaml.load(fh)

        yaml_content["channels"] = ["bioconda", "default"]

        with open(env_path, "w") as fh:
            fh.write(yaml_dump_to_string(yaml_content))

        module_lint = nf_core.modules.lint.ModuleLint(directory=self.nfcore_modules)
        module_lint.lint(print_results=False, module="bpipe/test")

        assert len(module_lint.failed) == 1, f"Linting failed with {[x.__dict__ for x in module_lint.failed]}"
        assert len(module_lint.passed) > 0
        assert len(module_lint.warned) >= 0
        assert module_lint.failed[0].lint_test == "environment_yml_valid"
