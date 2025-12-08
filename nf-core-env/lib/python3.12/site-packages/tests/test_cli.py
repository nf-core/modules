"""Tests covering the command-line code.

Most tests check the cli arguments are passed along and that some action is
taken.
"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from click.testing import CliRunner

import nf_core.__main__
import nf_core.utils


@mock.patch("nf_core.__main__.nf_core_cli")
def test_header(mock_cli):
    """Just try to execute the header function"""
    nf_core.__main__.run_nf_core()


@mock.patch("nf_core.__main__.nf_core_cli")
@mock.patch("nf_core.__main__.check_if_outdated", return_value=(True, None, "dummy_version"))
def test_header_outdated(mock_check_outdated, mock_nf_core_cli, capsys):
    """Check cli notifies the user when nf_core is outdated"""
    nf_core.__main__.run_nf_core()
    captured = capsys.readouterr()
    assert "There is a new version of nf-core/tools available! (dummy_version)" in captured.err


class TestCli(unittest.TestCase):
    """Class for testing the command line interface"""

    def setUp(self):
        self.runner = CliRunner()

    def assemble_params(self, params):
        """Assemble a dictionary of parameters into a list of arguments for the cli

        Note:
            if the value of a parameter is None, it will be considered a flag.
            Booleans were not used to avoid conflicting with the click.BOOL type.

        Args:
            params (dict): dict of parameters to assemble"""
        arg_list = []
        for key, value in params.items():
            if value is not None:
                arg_list += [f"--{key}", value]
            else:
                arg_list += [f"--{key}"]

        return arg_list

    def invoke_cli(self, cmd):
        """Invoke the commandline interface using a list of parameters

        Args:
            cmd (list): commandline to execute
        """
        return self.runner.invoke(nf_core.__main__.nf_core_cli, cmd)

    def test_cli_help(self):
        """Test the main launch function with --help"""
        result = self.invoke_cli(["--help"])
        assert result.exit_code == 0
        assert "Show the version and exit." in result.output

    def test_cli_bad_subcommand(self):
        """Test the main launch function with an unrecognised argument"""
        result = self.invoke_cli(["foo"])
        assert result.exit_code == 2

    def test_cli_verbose(self):
        """Test the main launch function with verbose flag"""
        result = self.invoke_cli(["-v"])
        # Checks that -v was considered valid
        assert "No such option: -v" not in nf_core.utils.strip_ansi_codes(result.output)

    @mock.patch("nf_core.pipelines.list.list_workflows", return_value="pipeline test list")
    def test_cli_list(self, mock_list_workflows):
        """Test nf-core pipelines are listed and cli parameters are passed on."""
        params = {
            "sort": "name",
            "json": None,
            "show-archived": None,
        }
        cmd = ["pipelines", "list"] + self.assemble_params(params) + ["kw1", "kw2"]
        result = self.invoke_cli(cmd)

        mock_list_workflows.assert_called_once_with(
            tuple(cmd[-2:]), params["sort"], "json" in params, "show-archived" in params
        )
        assert result.exit_code == 0
        assert "pipeline test list" in result.output

    @mock.patch("nf_core.pipelines.launch.Launch")
    def test_cli_launch(self, mock_launcher):
        """Test nf-core pipeline is launched and cli parameters are passed on."""
        mock_launcher.return_value.launch_pipeline.return_value = True

        temp_params_in = tempfile.NamedTemporaryFile()
        params = {
            "revision": "abcdef",
            "id": "idgui",
            "command-only": None,
            "params-out": "/path/params/out",
            "params-in": temp_params_in.name,
            "save-all": None,
            "show-hidden": None,
            "url": "builder_url",
        }
        cmd = ["pipelines", "launch"] + self.assemble_params(params) + ["pipeline_name"]
        result = self.invoke_cli(cmd)

        assert result.exit_code == 0

        mock_launcher.assert_called_once_with(
            cmd[-1],
            params["revision"],
            "command-only" in params,
            params["params-in"],
            params["params-out"],
            "save-all" in params,
            "show-hidden" in params,
            params["url"],
            params["id"],
        )

        mock_launcher.return_value.launch_pipeline.assert_called_once()

    @mock.patch("nf_core.pipelines.launch.Launch")
    def test_cli_launch_no_params_in(self, mock_launcher):
        """Test nf-core pipeline fails when params-in does not exist"""
        mock_launcher.return_value.launch_pipeline.return_value = True

        params = {
            "params-in": "/fake/path",
        }
        cmd = ["pipelines", "launch"] + self.assemble_params(params) + ["pipeline_name"]
        result = self.invoke_cli(cmd)

        assert result.exit_code == 2
        assert (
            f"Invalid value for '-p' / '--params-in': Path '{params['params-in']}' does not exist."
            in nf_core.utils.strip_ansi_codes(result.output)
        )

        mock_launcher.assert_not_called()

    @mock.patch("nf_core.pipelines.launch.Launch")
    def test_cli_launch_fail(self, mock_launcher):
        """Test nf-core pipeline fails with exit code 1 when pipeline fails."""
        mock_launcher.return_value.launch_pipeline.return_value = False
        cmd = ["pipelines", "launch", "pipeline_name"]
        result = self.invoke_cli(cmd)
        assert result.exit_code == 1

    @mock.patch("nf_core.pipelines.download.DownloadWorkflow")
    def test_cli_download(self, mock_dl):
        """Test nf-core pipeline is downloaded and cli parameters are passed on."""
        toplevel_params = {"hide-progress": None}
        params = {
            "revision": "abcdef",
            "outdir": "/path/outdir",
            "compress": "tar.gz",
            "force": None,
            "platform": None,
            "download-configuration": "yes",
            "tag": "3.14=testing",
            "container-system": "singularity",
            "container-library": "quay.io",
            "container-cache-utilisation": "copy",
            "container-cache-index": "/path/index.txt",
            "parallel-downloads": 2,
        }

        cmd = (
            self.assemble_params(toplevel_params)
            + ["pipelines", "download"]
            + self.assemble_params(params)
            + ["pipeline_name"]
        )
        result = self.invoke_cli(cmd)

        assert result.exit_code == 0

        mock_dl.assert_called_once_with(
            cmd[-1],
            (params["revision"],),
            params["outdir"],
            compress_type=params["compress"],
            force="force" in params,
            platform="platform" in params,
            download_configuration=params["download-configuration"],
            additional_tags=(params["tag"],),
            container_system=params["container-system"],
            container_library=(params["container-library"],),
            container_cache_utilisation=params["container-cache-utilisation"],
            container_cache_index=params["container-cache-index"],
            parallel=params["parallel-downloads"],
            hide_progress="hide-progress" in toplevel_params,
        )

        mock_dl.return_value.download_workflow.assert_called_once()

    @mock.patch("nf_core.pipelines.create.create.PipelineCreate")
    def test_create(self, mock_create):
        """Test nf-core pipeline is created and cli parameters are passed on."""
        params = {
            "name": "pipelinename",
            "description": "pipeline description",
            "author": "Kalle Anka",
            "outdir": "/path/outdir",
        }

        cmd = ["pipelines", "create"] + self.assemble_params(params)
        result = self.invoke_cli(cmd)

        assert result.exit_code == 0
        mock_create.assert_called_once_with(
            params["name"],
            params["description"],
            params["author"],
            force="force" in params,
            version="1.0.0dev",
            outdir=params["outdir"],
            template_config=None,
            organisation="nf-core",
        )
        mock_create.return_value.init_pipeline.assert_called_once()

    @mock.patch("nf_core.pipelines.create.create.PipelineCreate")
    def test_create_error(self, mock_create):
        """Test `nf-core pipelines create` run without providing all the arguments thorws an error."""
        params = {
            "name": "pipelinename",
        }

        cmd = ["pipelines", "create"] + self.assemble_params(params)
        result = self.invoke_cli(cmd)

        assert result.exit_code == 1
        assert "Partial arguments supplied." in result.output

    @mock.patch("nf_core.pipelines.create.PipelineCreateApp")
    def test_create_app(self, mock_create):
        """Test `nf-core pipelines create` runs an App."""
        cmd = ["pipelines", "create"]
        result = self.invoke_cli(cmd)

        assert result.return_value == (0 or None)
        assert "Launching interactive nf-core pipeline creation tool." in result.output

        mock_create.assert_called_once_with()
        mock_create.return_value.run.assert_called_once()

    @mock.patch("nf_core.utils.is_pipeline_directory")
    @mock.patch("nf_core.pipelines.lint.run_linting")
    def test_lint(self, mock_lint, mock_is_pipeline):
        """Test nf-core pipelines lint"""
        mock_lint_results = (mock.MagicMock, mock.MagicMock, mock.MagicMock)
        mock_lint_results[0].failed = []
        mock_lint_results[1].failed = []
        mock_lint_results[2].failed = []
        mock_lint.return_value = mock_lint_results

        temp_pipeline_dir = tempfile.NamedTemporaryFile()
        params = {
            "dir": temp_pipeline_dir.name,
            "release": None,
            "fix": "fix test",
            "key": "key test",
            "show-passed": None,
            "fail-ignored": None,
            "fail-warned": None,
            "markdown": "output_file.md",
            "json": "output_file.json",
        }

        cmd = ["pipelines", "lint"] + self.assemble_params(params)
        result = self.invoke_cli(cmd)

        assert result.exit_code == 0
        mock_lint.assert_called_once_with(
            params["dir"],
            "release" in params,
            (params["fix"],),
            (params["key"],),
            "show-passed" in params,
            "fail-ignored" in params,
            "fail-warned" in params,
            "test",
            params["markdown"],
            params["json"],
            "hide-progress" in params,
        )

    def test_lint_no_dir(self):
        """Test nf-core pipelines lint fails if --dir does not exist"""
        params = {
            "dir": "/bad/path",
        }

        cmd = ["pipelines", "lint"] + self.assemble_params(params)
        result = self.invoke_cli(cmd)

        assert result.exit_code == 2
        assert (
            f"Invalid value for '-d' / '--dir': Path '{params['dir']}' does not exist."
            in nf_core.utils.strip_ansi_codes(result.output)
        )

    @mock.patch("nf_core.utils.is_pipeline_directory")
    def test_lint_dir_is_not_pipeline(self, mock_is_pipeline):
        """Test nf-core pipelines lint logs an error if not called from a pipeline directory."""
        error_txt = "UserWarning has been raised"
        mock_is_pipeline.side_effect = UserWarning(error_txt)

        cmd = ["pipelines", "lint"]
        with self.assertLogs() as captured_logs:
            result = self.invoke_cli(cmd)

        assert result.exit_code == 1
        assert error_txt in captured_logs.output[-1]
        assert captured_logs.records[-1].levelname == "ERROR"

    @mock.patch("nf_core.utils.is_pipeline_directory")
    @mock.patch("nf_core.pipelines.lint.run_linting")
    def test_lint_log_assert_error(self, mock_lint, mock_is_pipeline):
        """Test nf-core pipelines lint logs assertion errors"""
        error_txt = "AssertionError has been raised"
        mock_lint.side_effect = AssertionError(error_txt)

        cmd = ["pipelines", "lint"]
        with self.assertLogs() as captured_logs:
            result = self.invoke_cli(cmd)

        assert result.exit_code == 1
        assert error_txt in captured_logs.output[-1]
        assert captured_logs.records[-1].levelname == "CRITICAL"

    @mock.patch("nf_core.utils.is_pipeline_directory")
    @mock.patch("nf_core.pipelines.lint.run_linting")
    def test_lint_log_user_warning(self, mock_lint, mock_is_pipeline):
        """Test nf-core pipelines lint logs assertion errors"""
        error_txt = "AssertionError has been raised"
        mock_lint.side_effect = UserWarning(error_txt)

        cmd = ["pipelines", "lint"]
        with self.assertLogs() as captured_logs:
            result = self.invoke_cli(cmd)

        assert result.exit_code == 1
        assert error_txt in captured_logs.output[-1]
        assert captured_logs.records[-1].levelname == "ERROR"

    @mock.patch("nf_core.pipelines.schema.PipelineSchema.get_schema_path")
    def test_schema_lint(self, mock_get_schema_path):
        """Test nf-core pipelines schema lint defaults to nextflow_schema.json"""
        cmd = ["pipelines", "schema", "lint"]
        with self.runner.isolated_filesystem():
            with open("nextflow_schema.json", "w") as f:
                f.write("{}")
            self.invoke_cli(cmd)
            mock_get_schema_path.assert_called_with(Path("nextflow_schema.json"))

    @mock.patch("nf_core.pipelines.schema.PipelineSchema.get_schema_path")
    def test_schema_lint_filename(self, mock_get_schema_path):
        """Test nf-core pipelines schema lint accepts a filename"""
        cmd = ["pipelines", "schema", "lint", "some_other_filename"]
        with self.runner.isolated_filesystem():
            with open("some_other_filename", "w") as f:
                f.write("{}")
            self.invoke_cli(cmd)
            mock_get_schema_path.assert_called_with(Path("some_other_filename"))

    @mock.patch("nf_core.pipelines.create_logo.create_logo")
    def test_create_logo(self, mock_create_logo):
        # Set up the mock to return a specific value

        cmd = ["pipelines", "create-logo", "test"]
        result = self.invoke_cli(cmd)

        mock_create_logo.assert_called_with("test", Path.cwd(), None, "light", 2300, "png", False)
        assert result.exit_code == 0
