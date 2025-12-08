#!/usr/bin/env python
"""nf-core: Helper tools for use with nf-core Nextflow pipelines."""

import logging
import os
import sys
from pathlib import Path

import rich
import rich.console
import rich.logging
import rich.traceback
import rich_click as click
import rich_click.rich_click as rc
from trogon import tui

from nf_core import __version__
from nf_core.commands_modules import (
    modules_bump_versions,
    modules_create,
    modules_info,
    modules_install,
    modules_lint,
    modules_list_local,
    modules_list_remote,
    modules_patch,
    modules_remove,
    modules_test,
    modules_update,
)
from nf_core.commands_pipelines import (
    pipelines_bump_version,
    pipelines_create,
    pipelines_create_logo,
    pipelines_create_params_file,
    pipelines_download,
    pipelines_launch,
    pipelines_lint,
    pipelines_list,
    pipelines_rocrate,
    pipelines_schema_build,
    pipelines_schema_docs,
    pipelines_schema_lint,
    pipelines_schema_validate,
    pipelines_sync,
)
from nf_core.commands_subworkflows import (
    subworkflows_create,
    subworkflows_info,
    subworkflows_install,
    subworkflows_lint,
    subworkflows_list_local,
    subworkflows_list_remote,
    subworkflows_remove,
    subworkflows_test,
    subworkflows_update,
)
from nf_core.commands_test_datasets import test_datasets_list_branches, test_datasets_list_remote, test_datasets_search
from nf_core.components.components_completion import autocomplete_modules, autocomplete_subworkflows
from nf_core.components.constants import NF_CORE_MODULES_REMOTE
from nf_core.pipelines.download.download import DownloadError
from nf_core.pipelines.list import autocomplete_pipelines
from nf_core.utils import check_if_outdated, nfcore_logo, rich_force_colors, setup_nfcore_dir

# Set up logging as the root logger
# Submodules should all traverse back to this
log = logging.getLogger()

# Set up .nfcore directory for storing files between sessions
setup_nfcore_dir()

# Set up nicer formatting of click cli help messages
rc.MAX_WIDTH = 100
rc.USE_RICH_MARKUP = True
rc.COMMANDS_BEFORE_OPTIONS = True


# Set up rich stderr console
stderr = rich.console.Console(stderr=True, force_terminal=rich_force_colors())
stdout = rich.console.Console(force_terminal=rich_force_colors())

# Set up the rich traceback
rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)


# Define exceptions for which no traceback should be printed,
# because they are actually preliminary, but intended program terminations.
# (Custom exceptions are cleaner than `sys.exit(1)`, which we used before)
def selective_traceback_hook(exctype, value, traceback):
    if exctype in {DownloadError, UserWarning, ValueError}:  # extend set as needed
        log.error(value)
    else:
        # print the colored traceback for all other exceptions with rich as usual
        stderr.print(rich.traceback.Traceback.from_exception(exctype, value, traceback))


sys.excepthook = selective_traceback_hook


# Define callback function to normalize the case of click arguments,
# which is used to make the module/subworkflow names, provided by the
# user on the cli, case insensitive.
def normalize_case(ctx, param, component_name):
    if component_name is not None:
        return component_name.casefold()


def run_nf_core():
    # print nf-core header if environment variable is not set
    if os.environ.get("_NF_CORE_COMPLETE") is None:
        # Print nf-core header
        stderr.print("\n")
        for line in nfcore_logo:
            stderr.print(line, highlight=False)
        stderr.print(
            f"\n[grey39]    nf-core/tools version {__version__} - [link=https://nf-co.re]https://nf-co.re[/]",
            highlight=False,
        )
        try:
            is_outdated, _, remote_vers = check_if_outdated()
            if is_outdated:
                stderr.print(
                    f"[bold bright_yellow]    There is a new version of nf-core/tools available! ({remote_vers})",
                    highlight=False,
                )
        except Exception as e:
            log.debug(f"Could not check latest version: {e}")
        stderr.print("\n")
    # Launch the click cli
    nf_core_cli(auto_envvar_prefix="NFCORE")


@tui(command="interface", help="Launch the nf-core interface")
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Print verbose output to the console.",
)
@click.option("--hide-progress", is_flag=True, default=False, help="Don't show progress bars.")
@click.option("-l", "--log-file", help="Save a verbose log to a file.", metavar="<filename>")
@click.pass_context
@click.rich_config(
    {
        "theme": "default-nu",
        "options_table_column_types": ["opt_long", "opt_short", "help"],
    }
)
def nf_core_cli(ctx, verbose, hide_progress, log_file):
    """
    nf-core/tools provides a set of helper tools for use with nf-core Nextflow pipelines.

    It is designed for both end-users running pipelines and also developers creating new pipelines.
    """
    # Set the base logger to output DEBUG
    log.setLevel(logging.DEBUG)

    # Set up logs to the console
    log.addHandler(
        rich.logging.RichHandler(
            level=logging.DEBUG if verbose else logging.INFO,
            console=rich.console.Console(stderr=True, force_terminal=rich_force_colors()),
            show_time=False,
            show_path=verbose,  # True if verbose, false otherwise
            markup=True,
        )
    )

    # don't show rich debug logging in verbose mode
    rich_logger = logging.getLogger("rich")
    rich_logger.setLevel(logging.INFO)

    # Set up logs to a file if we asked for one
    if log_file:
        log_fh = logging.FileHandler(log_file, encoding="utf-8")
        log_fh.setLevel(logging.DEBUG)
        log_fh.setFormatter(logging.Formatter("[%(asctime)s] %(name)-20s [%(levelname)-7s]  %(message)s"))
        log.addHandler(log_fh)

    ctx.obj = {
        "verbose": verbose,
        "hide_progress": hide_progress or verbose,  # Always hide progress bar with verbose logging
    }


# nf-core pipelines subcommands
@nf_core_cli.group(aliases=["p", "pipeline"])
@click.command_panel("For users", commands=["download", "create-params-file", "launch", "list"])
@click.command_panel(
    "For developers", commands=["bump-version", "create", "create-logo", "lint", "rocrate", "schema", "sync"]
)
@click.pass_context
def pipelines(ctx):
    """
    Commands to manage nf-core pipelines.
    """
    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)


# nf-core pipelines create
@pipelines.command("create")
@click.pass_context
@click.option(
    "-n",
    "--name",
    type=str,
    help="The name of your new pipeline",
)
@click.option("-d", "--description", type=str, help="A short description of your pipeline")
@click.option("-a", "--author", type=str, help="Name of the main author(s)")
@click.option("--version", type=str, default="1.0.0dev", help="The initial version number to use")
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite output directory if it already exists")
@click.option("-o", "--outdir", help="Output directory for new pipeline (default: pipeline name)")
@click.option("-t", "--template-yaml", help="Pass a YAML file to customize the template")
@click.option(
    "--organisation",
    type=str,
    default="nf-core",
    help="The name of the GitHub organisation where the pipeline will be hosted (default: nf-core)",
)
def command_pipelines_create(ctx, name, description, author, version, force, outdir, template_yaml, organisation):
    """
    Create a new pipeline using the nf-core template.
    """
    pipelines_create(ctx, name, description, author, version, force, outdir, template_yaml, organisation)


# nf-core pipelines lint
@pipelines.command("lint")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory [dim]\[default: current working directory][/]",
)
@click.option(
    "--release",
    is_flag=True,
    default=Path(os.environ.get("GITHUB_REF", "").strip(" '\"")).parent.name in ["master", "main"]
    and os.environ.get("GITHUB_REPOSITORY", "").startswith("nf-core/")
    and not os.environ.get("GITHUB_REPOSITORY", "") == "nf-core/tools",
    help="Execute additional checks for release-ready workflows.",
)
@click.option(
    "-f",
    "--fix",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Attempt to automatically fix specified lint test",
)
@click.option(
    "-k",
    "--key",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Run only these lint tests",
)
@click.option("-p", "--show-passed", is_flag=True, help="Show passing tests on the command line")
@click.option("-i", "--fail-ignored", is_flag=True, help="Convert ignored tests to failures")
@click.option("-w", "--fail-warned", is_flag=True, help="Convert warn tests to failures")
@click.option(
    "--markdown",
    type=str,
    metavar="<filename>",
    help="File to write linting results to (Markdown)",
)
@click.option(
    "--json",
    type=str,
    metavar="<filename>",
    help="File to write linting results to (JSON)",
)
@click.option(
    "--sort-by",
    type=click.Choice(["module", "test"]),
    default="test",
    help="Sort lint output by module or test name.",
    show_default=True,
)
@click.pass_context
def command_pipelines_lint(
    ctx,
    directory,
    release,
    fix,
    key,
    show_passed,
    fail_ignored,
    fail_warned,
    markdown,
    json,
    sort_by,
):
    """
    Check pipeline code against nf-core guidelines.
    """
    pipelines_lint(ctx, directory, release, fix, key, show_passed, fail_ignored, fail_warned, markdown, json, sort_by)


# nf-core pipelines download
@pipelines.command("download")
@click.argument(
    "pipeline",
    required=False,
    metavar="<pipeline name>",
    shell_complete=autocomplete_pipelines,
)
@click.option(
    "-r",
    "--revision",
    multiple=True,
    help="Pipeline release to download. Multiple invocations are possible, e.g. `-r 1.1 -r 1.2`",
)
@click.option("-o", "--outdir", type=str, help="Output directory")
@click.option(
    "-x",
    "--compress",
    type=click.Choice(["tar.gz", "tar.bz2", "zip", "none"]),
    help="Archive compression type",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-p",
    "--platform",
    is_flag=True,
    default=False,
    help="Download for Seqera Platform (formerly Nextflow Tower)",
)
@click.option(
    "-c",
    "--download-configuration",
    type=click.Choice(["yes", "no"]),
    default="no",
    help="Include configuration profiles in download. Not available with `--platform`",
)
@click.option(
    "--tag",
    multiple=True,
    help="Add custom alias tags to `--platform` downloads. For example, `--tag \"3.10=validated\"` adds the custom 'validated' tag to the 3.10 release.",
)
# -c changed to -s for consistency with other --container arguments, where it is always the first letter of the last word.
# Also -c might be used instead of -d for config in a later release, but reusing params for different options in two subsequent releases might be too error-prone.
@click.option(
    "-s",
    "--container-system",
    type=click.Choice(["none", "singularity", "docker"]),
    help="Download container images of required software.",
)
@click.option(
    "-l",
    "--container-library",
    multiple=True,
    help="Container registry/library or mirror to pull images from. Not available for Docker containers.",
)
@click.option(
    "-u",
    "--container-cache-utilisation",
    type=click.Choice(["amend", "copy", "remote"]),
    help="Utilise a `singularity.cacheDir` in the download process, if applicable. Not available for Docker containers.",
)
@click.option(
    "-i",
    "--container-cache-index",
    type=str,
    help="List of images already available in a remote `singularity.cacheDir`. Not available for Docker containers.",
)
@click.option(
    "-d",
    "--parallel-downloads",
    type=int,
    default=4,
    help="Number of allowed parallel tasks",
)
@click.pass_context
def command_pipelines_download(
    ctx,
    pipeline,
    revision,
    outdir,
    compress,
    force,
    platform,
    download_configuration,
    tag,
    container_system,
    container_library,
    container_cache_utilisation,
    container_cache_index,
    parallel_downloads,
):
    """
    Download a pipeline, nf-core/configs and pipeline singularity images.
    """
    pipelines_download(
        ctx,
        pipeline,
        revision,
        outdir,
        compress,
        force,
        platform,
        download_configuration,
        tag,
        container_system,
        container_library,
        container_cache_utilisation,
        container_cache_index,
        parallel_downloads,
    )


# nf-core pipelines create-params-file
@pipelines.command("create-params-file")
@click.argument(
    "pipeline",
    required=False,
    metavar="<pipeline name>",
    shell_complete=autocomplete_pipelines,
)
@click.option("-r", "--revision", help="Release/branch/SHA of the pipeline (if remote)")
@click.option(
    "-o",
    "--output",
    type=str,
    default="nf-params.yml",
    metavar="<filename>",
    help="Output filename. Defaults to `nf-params.yml`.",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-x",
    "--show-hidden",
    is_flag=True,
    default=False,
    help="Show hidden params which don't normally need changing",
)
@click.pass_context
def command_pipelines_create_params_file(ctx, pipeline, revision, output, force, show_hidden):
    """
    Build a parameter file for a pipeline.
    """
    pipelines_create_params_file(ctx, pipeline, revision, output, force, show_hidden)


# nf-core pipelines launch
@pipelines.command("launch")
@click.argument(
    "pipeline",
    required=False,
    metavar="<pipeline name>",
    shell_complete=autocomplete_pipelines,
)
@click.option("-r", "--revision", help="Release/branch/SHA of the project to run (if remote)")
@click.option("-i", "--id", help="ID for web-gui launch parameter set")
@click.option(
    "-c",
    "--command-only",
    is_flag=True,
    default=False,
    help="Create Nextflow command with params (no params file)",
)
@click.option(
    "-o",
    "--params-out",
    type=click.Path(),
    default=os.path.join(os.getcwd(), "nf-params.json"),
    help="Path to save run parameters file",
)
@click.option(
    "-p",
    "--params-in",
    type=click.Path(exists=True),
    help="Set of input run params to use from a previous run",
)
@click.option(
    "-a",
    "--save-all",
    is_flag=True,
    default=False,
    help="Save all parameters, even if unchanged from default",
)
@click.option(
    "-x",
    "--show-hidden",
    is_flag=True,
    default=False,
    help="Show hidden params which don't normally need changing",
)
@click.option(
    "-u",
    "--url",
    type=str,
    default="https://nf-co.re/launch",
    help="Customise the builder URL (for development work)",
)
@click.pass_context
def command_pipelines_launch(
    ctx,
    pipeline,
    id,
    revision,
    command_only,
    params_in,
    params_out,
    save_all,
    show_hidden,
    url,
):
    """
    Launch a pipeline using a web GUI or command line prompts.
    """
    pipelines_launch(ctx, pipeline, id, revision, command_only, params_in, params_out, save_all, show_hidden, url)


# nf-core pipelines list
@pipelines.command("list")
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option(
    "-s",
    "--sort",
    type=click.Choice(["release", "pulled", "name", "stars"]),
    default="release",
    help="How to sort listed pipelines",
)
@click.option("--json", is_flag=True, default=False, help="Print full output as JSON")
@click.option("--show-archived", is_flag=True, default=False, help="Print archived workflows")
@click.pass_context
def command_pipelines_list(ctx, keywords, sort, json, show_archived):
    """
    List available nf-core pipelines with local info.
    """
    pipelines_list(ctx, keywords, sort, json, show_archived)


# nf-core pipelines rocrate
@pipelines.command("rocrate")
@click.argument(
    "pipeline_dir",
    type=click.Path(exists=True),
    default=Path.cwd(),
    required=True,
    metavar="<pipeline directory>",
)
@click.option(
    "-j",
    "--json_path",
    default=Path.cwd(),
    type=str,
    help="Path to save RO Crate metadata json file to",
)
@click.option("-z", "--zip_path", type=str, help="Path to save RO Crate zip file to")
@click.option(
    "-pv",
    "--pipeline_version",
    type=str,
    help="Version of pipeline to use for RO Crate",
    default="",
)
@click.pass_context
def rocrate(
    ctx,
    pipeline_dir: str,
    json_path: str,
    zip_path: str,
    pipeline_version: str,
):
    """
    Make an Research Object Crate
    """
    pipelines_rocrate(ctx, pipeline_dir, json_path, zip_path, pipeline_version)


# nf-core pipelines sync
@pipelines.command("sync")
@click.pass_context
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-b",
    "--from-branch",
    type=str,
    help="The git branch to use to fetch workflow variables.",
)
@click.option(
    "-p",
    "--pull-request",
    is_flag=True,
    default=False,
    help="Make a GitHub pull-request with the changes.",
)
@click.option(
    "--force_pr",
    is_flag=True,
    default=False,
    help="Force the creation of a pull-request, even if there are no changes.",
)
@click.option("-g", "--github-repository", type=str, help="GitHub PR: target repository.")
@click.option("-u", "--username", type=str, help="GitHub PR: auth username.")
@click.option("-t", "--template-yaml", help="Pass a YAML file to customize the template")
@click.option("-b", "--blog-post", type=str, help="Link to the blog post")
def command_pipelines_sync(
    ctx, directory, from_branch, pull_request, github_repository, username, template_yaml, force_pr, blog_post
):
    """
    Sync a pipeline [cyan i]TEMPLATE[/] branch with the nf-core template.
    """
    pipelines_sync(
        ctx, directory, from_branch, pull_request, github_repository, username, template_yaml, force_pr, blog_post
    )


# nf-core pipelines bump-version
@pipelines.command("bump-version")
@click.pass_context
@click.argument("new_version", required=True, metavar="<new version>")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-n",
    "--nextflow",
    is_flag=True,
    default=False,
    help="Bump required nextflow version instead of pipeline version",
)
def command_pipelines_bump_version(ctx, new_version, directory, nextflow):
    """
    Update nf-core pipeline version number with `nf-core pipelines bump-version`.
    """
    pipelines_bump_version(ctx, new_version, directory, nextflow)


# nf-core pipelines create-logo
@pipelines.command("create-logo")
@click.argument("logo-text", metavar="<logo_text>")
@click.option("-d", "--dir", "directory", type=click.Path(), default=".", help="Directory to save the logo in.")
@click.option(
    "-n",
    "--name",
    type=str,
    help="Name of the output file (with or without '.png' suffix).",
)
@click.option(
    "--theme",
    type=click.Choice(["light", "dark"]),
    default="light",
    help="Theme for the logo.",
    show_default=True,
)
@click.option(
    "--width",
    type=int,
    default=2300,
    help="Width of the logo in pixels.",
    show_default=True,
)
@click.option(
    "--format",
    type=click.Choice(["png", "svg"]),
    default="png",
    help="Image format of the logo, either PNG or SVG.",
    show_default=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite any files if they already exist",
)
def command_pipelines_create_logo(logo_text, directory, name, theme, width, format, force):
    """
    Generate a logo with the nf-core logo template.
    """
    pipelines_create_logo(logo_text, directory, name, theme, width, format, force)


# nf-core pipelines schema subcommands
@pipelines.group("schema")
def pipeline_schema():
    """
    Suite of tools for developers to manage pipeline schema.

    All nf-core pipelines should have a nextflow_schema.json file in their
    root directory that describes the different pipeline parameters.
    """
    pass


# nf-core pipelines schema validate
@pipeline_schema.command("validate")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.argument(
    "pipeline",
    required=False,
    metavar="<pipeline name>",
    shell_complete=autocomplete_pipelines,
)
@click.argument("params", type=click.Path(exists=True), required=True, metavar="<JSON params file>")
def command_pipelines_schema_validate(directory, pipeline, params):
    """
    Validate a set of parameters against a pipeline schema.
    """
    if Path(directory, pipeline).exists():
        # this is a local pipeline
        pipeline = Path(directory, pipeline)

    pipelines_schema_validate(pipeline, params)


# nf-core pipelines schema build
@pipeline_schema.command("build")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "--no-prompts",
    is_flag=True,
    help="Do not confirm changes, just update parameters and exit",
)
@click.option(
    "--web-only",
    is_flag=True,
    help="Skip building using Nextflow config, just launch the web tool",
)
@click.option(
    "--url",
    type=str,
    default="https://nf-co.re/pipeline_schema_builder",
    help="Customise the builder URL (for development work)",
)
def command_pipelines_schema_build(directory, no_prompts, web_only, url):
    """
    Interactively build a pipeline schema from Nextflow params.
    """
    pipelines_schema_build(directory, no_prompts, web_only, url)


# nf-core pipelines schema lint
@pipeline_schema.command("lint")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.argument(
    "schema_file",
    type=click.Path(exists=True),
    default="nextflow_schema.json",
    metavar="<pipeline schema>",
)
def command_pipelines_schema_lint(directory, schema_file):
    """
    Check that a given pipeline schema is valid.
    """
    pipelines_schema_lint(Path(directory, schema_file))


# nf-core pipelines schema docs
@pipeline_schema.command("docs")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.argument(
    "schema_file",
    type=click.Path(exists=True),
    default="nextflow_schema.json",
    required=False,
    metavar="<pipeline schema>",
)
@click.option(
    "-o",
    "--output",
    type=str,
    metavar="<filename>",
    help="Output filename. Defaults to standard out.",
)
@click.option(
    "-x",
    "--format",
    type=click.Choice(["markdown", "html"]),
    default="markdown",
    help="Format to output docs in.",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-c",
    "--columns",
    type=str,
    metavar="<columns_list>",
    help="CSV list of columns to include in the parameter tables (parameter,description,type,default,required,hidden)",
    default="parameter,description,type,default,required,hidden",
)
def command_pipelines_schema_docs(directory, schema_file, output, format, force, columns):
    """
    Outputs parameter documentation for a pipeline schema.
    """
    pipelines_schema_docs(Path(directory, schema_file), output, format, force, columns)


# nf-core modules subcommands
@nf_core_cli.group(aliases=["m", "module"])
@click.option(
    "-g",
    "--git-remote",
    type=str,
    default=NF_CORE_MODULES_REMOTE,
    help="Remote git repo to fetch files from",
)
@click.option(
    "-b",
    "--branch",
    type=str,
    default=None,
    help="Branch of git repository hosting modules.",
)
@click.option(
    "-N",
    "--no-pull",
    is_flag=True,
    default=False,
    help="Do not pull in latest changes to local clone of modules repository.",
)
@click.command_panel("For pipeline development", commands=["list", "info", "install", "update", "remove", "patch"])
@click.command_panel("For module development", commands=["create", "lint", "test", "bump-versions"])
@click.pass_context
def modules(ctx, git_remote, branch, no_pull):
    """
    Commands to manage Nextflow DSL2 modules (tool wrappers).
    """
    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)

    # Place the arguments in a context object
    ctx.obj["modules_repo_url"] = git_remote
    ctx.obj["modules_repo_branch"] = branch
    ctx.obj["modules_repo_no_pull"] = no_pull


# nf-core modules list subcommands
@modules.group("list")
@click.pass_context
def modules_list(ctx):
    """
    List modules in a local pipeline or remote repository.
    """
    pass


# nf-core modules list remote
@modules_list.command("remote")
@click.pass_context
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option("-j", "--json", is_flag=True, help="Print as JSON to stdout")
def command_modules_list_remote(ctx, keywords, json):
    """
    List modules in a remote GitHub repo [dim i](e.g [link=https://github.com/nf-core/modules]nf-core/modules[/])[/].
    """
    modules_list_remote(ctx, keywords, json)


# nf-core modules list local
@modules_list.command("local")
@click.pass_context
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option("-j", "--json", is_flag=True, help="Print as JSON to stdout")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: Current working directory][/]",
)
def command_modules_list_local(ctx, keywords, json, directory):  # pylint: disable=redefined-builtin
    """
    List modules installed locally in a pipeline
    """
    modules_list_local(ctx, keywords, json, directory)


# nf-core modules install
@modules.command("install")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-p",
    "--prompt",
    is_flag=True,
    default=False,
    help="Prompt for the version of the module",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Force reinstallation of module if it already exists",
)
@click.option("-s", "--sha", type=str, metavar="<commit sha>", help="Install module at commit SHA")
def command_modules_install(ctx, tool, directory, prompt, force, sha):
    """
    Install DSL2 modules within a pipeline.
    """
    modules_install(ctx, tool, directory, prompt, force, sha)


# nf-core modules update
@modules.command("update")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Force update of module")
@click.option(
    "-p",
    "--prompt",
    is_flag=True,
    default=False,
    help="Prompt for the version of the module",
)
@click.option(
    "--limit-output",
    "limit_output",
    is_flag=True,
    default=False,
    help="Limit output to only the difference in main.nf",
)
@click.option("-s", "--sha", type=str, metavar="<commit sha>", help="Install module at commit SHA")
@click.option(
    "-a",
    "--all",
    "install_all",
    is_flag=True,
    default=False,
    help="Update all modules installed in pipeline",
)
@click.option(
    "-x/-y",
    "--preview/--no-preview",
    is_flag=True,
    default=None,
    help="Preview / no preview of changes before applying",
)
@click.option(
    "-D",
    "--save-diff",
    type=str,
    metavar="<filename>",
    default=None,
    help="Save diffs to a file instead of updating in place",
)
@click.option(
    "-u",
    "--update-deps",
    is_flag=True,
    default=False,
    help="Automatically update all linked modules and subworkflows without asking for confirmation",
)
def command_modules_update(
    ctx,
    tool,
    directory,
    force,
    prompt,
    sha,
    install_all,
    preview,
    save_diff,
    update_deps,
    limit_output,
):
    """
    Update DSL2 modules within a pipeline.
    """
    modules_update(ctx, tool, directory, force, prompt, sha, install_all, preview, save_diff, update_deps, limit_output)


# nf-core modules patch
@modules.command("patch")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option("-r", "--remove", is_flag=True, default=False, help="Remove an existent patch file and regenerate it.")
def command_modules_patch(ctx, tool, directory, remove):
    """
    Create a patch file for minor changes in a module
    """
    modules_patch(ctx, tool, directory, remove)


# nf-core modules remove
@modules.command("remove")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
def command_modules_remove(ctx, directory, tool):
    """
    Remove a module from a pipeline.
    """
    modules_remove(ctx, directory, tool)


# nf-core modules create
@modules.command("create")
@click.pass_context
@click.argument("tool", type=str, required=False, metavar="<tool> or <tool/subtool>")
@click.option("-d", "--dir", "directory", type=click.Path(exists=True), default=".", metavar="<directory>")
@click.option(
    "-a",
    "--author",
    type=str,
    metavar="<author>",
    help="Module author's GitHub username prefixed with '@'",
)
@click.option(
    "-l",
    "--label",
    type=str,
    metavar="<process label>",
    help="Standard resource label for process",
)
@click.option(
    "-m",
    "--meta",
    is_flag=True,
    default=False,
    help="Use Groovy meta map for sample information",
)
@click.option(
    "-n",
    "--no-meta",
    is_flag=True,
    default=False,
    help="Don't use meta map for sample information",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite any files if they already exist",
)
@click.option(
    "-c",
    "--conda-name",
    type=str,
    default=None,
    help="Name of the conda package to use",
)
@click.option(
    "-p",
    "--conda-package-version",
    type=str,
    default=None,
    help="Version of conda package to use",
)
@click.option(
    "-i",
    "--empty-template",
    is_flag=True,
    default=False,
    help="Create a module from the template without TODOs or examples",
)
@click.option(
    "--migrate-pytest",
    is_flag=True,
    default=False,
    help="Migrate a module with pytest tests to nf-test",
)
def command_modules_create(
    ctx,
    tool,
    directory,
    author,
    label,
    meta,
    no_meta,
    force,
    conda_name,
    conda_package_version,
    empty_template,
    migrate_pytest,
):
    """
    Create a new DSL2 module from the nf-core template.
    """
    modules_create(
        ctx,
        tool,
        directory,
        author,
        label,
        meta,
        no_meta,
        force,
        conda_name,
        conda_package_version,
        empty_template,
        migrate_pytest,
    )


# nf-core modules test
@modules.command("test")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Print verbose output to the console. Sets `--debug` inside the nf-test command.",
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    metavar="<nf-core/modules directory>",
)
@click.option(
    "-p",
    "--no-prompts",
    is_flag=True,
    default=False,
    help="Use defaults without prompting",
)
@click.option("-u", "--update", is_flag=True, default=False, help="Update existing snapshots")
@click.option(
    "-o",
    "--once",
    is_flag=True,
    default=False,
    help="Run tests only once. Don't check snapshot stability",
)
@click.option(
    "--profile",
    type=click.Choice(["docker", "singularity", "conda"]),
    default=None,
    help="Run tests with a specific profile",
)
@click.option(
    "--migrate-pytest",
    is_flag=True,
    default=False,
    help="Migrate a module with pytest tests to nf-test",
)
def command_modules_test(ctx, tool, directory, no_prompts, update, once, profile, migrate_pytest, verbose):
    """
    Run nf-test for a module.
    """
    if verbose:
        ctx.obj["verbose"] = verbose
    modules_test(ctx, tool, directory, no_prompts, update, once, profile, migrate_pytest)


# nf-core modules lint
@modules.command("lint")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    metavar="<pipeline/modules directory>",
)
@click.option(
    "-r",
    "--registry",
    type=str,
    metavar="<registry>",
    default=None,
    help="Registry to use for containers. If not specified it will use docker.registry value in the nextflow.config file",
)
@click.option(
    "-k",
    "--key",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Run only these lint tests",
)
@click.option("-a", "--all", is_flag=True, help="Run on all modules")
@click.option("-w", "--fail-warned", is_flag=True, help="Convert warn tests to failures")
@click.option("--local", is_flag=True, help="Run additional lint tests for local modules")
@click.option("--passed", is_flag=True, help="Show passed tests")
@click.option(
    "--sort-by",
    type=click.Choice(["module", "test"]),
    default="test",
    help="Sort lint output by module or test name.",
    show_default=True,
)
@click.option(
    "--fix-version",
    is_flag=True,
    help="Fix the module version if a newer version is available",
)
@click.option("--fix", is_flag=True, help="Fix all linting tests if possible.")
def command_modules_lint(
    ctx, tool, directory, registry, key, all, fail_warned, local, passed, sort_by, fix_version, fix
):
    """
    Lint one or more modules in a directory.
    """
    modules_lint(ctx, tool, directory, registry, key, all, fail_warned, local, passed, sort_by, fix_version, fix)


# nf-core modules info
@modules.command("info")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: Current working directory][/]",
)
def command_modules_info(ctx, tool, directory):
    """
    Show developer usage information about a given module.
    """
    modules_info(ctx, tool, directory)


# nf-core modules bump-versions
@modules.command("bump-versions")
@click.pass_context
@click.argument(
    "tool",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="<tool> or <tool/subtool>. Module to bump versions for. If <tool> is provided and <tool/subtool> exists, all subtools will be bumped.",
    shell_complete=autocomplete_modules,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    metavar="<nf-core/modules directory>",
)
@click.option("-a", "--all", is_flag=True, help="Run on all modules")
@click.option("-s", "--show-all", is_flag=True, help="Show up-to-date modules in results too")
@click.option("-r", "--dry-run", is_flag=True, help="Dry run the command")
def command_modules_bump_versions(ctx, tool, directory, all, show_all, dry_run):
    """
    Bump versions for one or more modules in a clone of
    the nf-core/modules repo.
    """
    modules_bump_versions(ctx, tool, directory, all, show_all, dry_run)


# nf-core subworkflows click command
@nf_core_cli.group(aliases=["s", "swf", "subworkflow"])
@click.option(
    "-g",
    "--git-remote",
    type=str,
    default=NF_CORE_MODULES_REMOTE,
    help="Remote git repo to fetch files from",
)
@click.option(
    "-b",
    "--branch",
    type=str,
    default=None,
    help="Branch of git repository hosting modules.",
)
@click.option(
    "-N",
    "--no-pull",
    is_flag=True,
    default=False,
    help="Do not pull in latest changes to local clone of modules repository.",
)
@click.command_panel("For pipeline development", commands=["list", "info", "install", "update", "remove", "patch"])
@click.command_panel("For module development", commands=["create", "lint", "test", "bump-versions"])
@click.pass_context
def subworkflows(ctx, git_remote, branch, no_pull):
    """
    Commands to manage Nextflow DSL2 subworkflows (tool wrappers).
    """
    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)

    # Place the arguments in a context object
    ctx.obj["modules_repo_url"] = git_remote
    ctx.obj["modules_repo_branch"] = branch
    ctx.obj["modules_repo_no_pull"] = no_pull


# nf-core subworkflows create
@subworkflows.command("create")
@click.pass_context
@click.argument("subworkflow", type=str, required=False, metavar="subworkflow name")
@click.option("-d", "--dir", "directory", type=click.Path(exists=True), default=".", metavar="<directory>")
@click.option(
    "-a",
    "--author",
    type=str,
    metavar="<author>",
    help="Module author's GitHub username prefixed with '@'",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite any files if they already exist",
)
@click.option(
    "--migrate-pytest",
    is_flag=True,
    default=False,
    help="Migrate a module with pytest tests to nf-test",
)
def command_subworkflows_create(ctx, subworkflow, directory, author, force, migrate_pytest):
    """
    Create a new subworkflow from the nf-core template.
    """
    subworkflows_create(ctx, subworkflow, directory, author, force, migrate_pytest)


# nf-core subworkflows test
@subworkflows.command("test")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    metavar="<nf-core/modules directory>",
)
@click.option(
    "-p",
    "--no-prompts",
    is_flag=True,
    default=False,
    help="Use defaults without prompting",
)
@click.option("-u", "--update", is_flag=True, default=False, help="Update existing snapshots")
@click.option(
    "-o",
    "--once",
    is_flag=True,
    default=False,
    help="Run tests only once. Don't check snapshot stability",
)
@click.option(
    "--profile",
    type=click.Choice(["docker", "singularity", "conda"]),
    default=None,
    help="Run tests with a specific profile",
)
@click.option(
    "--migrate-pytest",
    is_flag=True,
    default=False,
    help="Migrate a subworkflow with pytest tests to nf-test",
)
def command_subworkflows_test(ctx, subworkflow, directory, no_prompts, update, once, profile, migrate_pytest):
    """
    Run nf-test for a subworkflow.
    """
    subworkflows_test(ctx, subworkflow, directory, no_prompts, update, once, profile, migrate_pytest)


# nf-core subworkflows list subcommands
@subworkflows.group("list")
@click.pass_context
def subworkflows_list(ctx):
    """
    List subworkflows in a local pipeline or remote repository.
    """
    pass


# nf-core subworkflows list remote
@subworkflows_list.command("remote")
@click.pass_context
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option("-j", "--json", is_flag=True, help="Print as JSON to stdout")
def command_subworkflows_list_remote(ctx, keywords, json):
    """
    List subworkflows in a remote GitHub repo [dim i](e.g [link=https://github.com/nf-core/modules]nf-core/modules[/])[/].
    """
    subworkflows_list_remote(ctx, keywords, json)


# nf-core subworkflows list local
@subworkflows_list.command("local")
@click.pass_context
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option("-j", "--json", is_flag=True, help="Print as JSON to stdout")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: Current working directory][/]",
)
def command_subworkflows_list_local(ctx, keywords, json, directory):  # pylint: disable=redefined-builtin
    """
    List subworkflows installed locally in a pipeline
    """
    subworkflows_list_local(ctx, keywords, json, directory)


# nf-core subworkflows lint
@subworkflows.command("lint")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    metavar="<pipeline/modules directory>",
)
@click.option(
    "-r",
    "--registry",
    type=str,
    metavar="<registry>",
    default=None,
    help="Registry to use for containers. If not specified it will use docker.registry value in the nextflow.config file",
)
@click.option(
    "-k",
    "--key",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Run only these lint tests",
)
@click.option("-a", "--all", is_flag=True, help="Run on all subworkflows")
@click.option("-w", "--fail-warned", is_flag=True, help="Convert warn tests to failures")
@click.option("--local", is_flag=True, help="Run additional lint tests for local subworkflows")
@click.option("--passed", is_flag=True, help="Show passed tests")
@click.option(
    "--sort-by",
    type=click.Choice(["subworkflow", "test"]),
    default="test",
    help="Sort lint output by subworkflow or test name.",
    show_default=True,
)
@click.option("--fix", is_flag=True, help="Fix all linting tests if possible.")
def command_subworkflows_lint(
    ctx, subworkflow, directory, registry, key, all, fail_warned, local, passed, sort_by, fix
):
    """
    Lint one or more subworkflows in a directory.
    """
    subworkflows_lint(ctx, subworkflow, directory, registry, key, all, fail_warned, local, passed, sort_by, fix)


# nf-core subworkflows info
@subworkflows.command("info")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: Current working directory][/]",
)
def command_subworkflows_info(ctx, subworkflow, directory):
    """
    Show developer usage information about a given subworkflow.
    """
    subworkflows_info(ctx, subworkflow, directory)


# nf-core subworkflows install
@subworkflows.command("install")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-p",
    "--prompt",
    is_flag=True,
    default=False,
    help="Prompt for the version of the subworkflow",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Force reinstallation of subworkflow if it already exists",
)
@click.option(
    "-s",
    "--sha",
    type=str,
    metavar="<commit sha>",
    help="Install subworkflow at commit SHA",
)
def command_subworkflows_install(ctx, subworkflow, directory, prompt, force, sha):
    """
    Install DSL2 subworkflow within a pipeline.
    """
    subworkflows_install(ctx, subworkflow, directory, prompt, force, sha)


# nf-core subworkflows patch
@subworkflows.command("patch")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option("-r", "--remove", is_flag=True, default=False, help="Remove an existent patch file and regenerate it.")
def subworkflows_patch(ctx, subworkflow, dir, remove):
    """
    Create a patch file for minor changes in a subworkflow

    Checks if a subworkflow has been modified locally and creates a patch file
    describing how the module has changed from the remote version
    """
    from nf_core.subworkflows import SubworkflowPatch

    try:
        subworkflow_patch = SubworkflowPatch(
            dir,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )
        if remove:
            subworkflow_patch.remove(subworkflow)
        else:
            subworkflow_patch.patch(subworkflow)
    except (UserWarning, LookupError) as e:
        log.error(e)
        sys.exit(1)


# nf-core subworkflows remove
@subworkflows.command("remove")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
def command_subworkflows_remove(ctx, directory, subworkflow):
    """
    Remove a subworkflow from a pipeline.
    """
    subworkflows_remove(ctx, directory, subworkflow)


# nf-core subworkflows update
@subworkflows.command("update")
@click.pass_context
@click.argument(
    "subworkflow",
    type=str,
    callback=normalize_case,
    required=False,
    metavar="subworkflow name",
    shell_complete=autocomplete_subworkflows,
)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Force update of subworkflow")
@click.option(
    "-p",
    "--prompt",
    is_flag=True,
    default=False,
    help="Prompt for the version of the subworkflow",
)
@click.option(
    "-s",
    "--sha",
    type=str,
    metavar="<commit sha>",
    help="Install subworkflow at commit SHA",
)
@click.option(
    "-l",
    "--limit-output",
    "limit_output",
    is_flag=True,
    default=False,
    help="Limit output to only the difference in main.nf",
)
@click.option(
    "-a",
    "--all",
    "install_all",
    is_flag=True,
    default=False,
    help="Update all subworkflow installed in pipeline",
)
@click.option(
    "-x/-y",
    "--preview/--no-preview",
    is_flag=True,
    default=None,
    help="Preview / no preview of changes before applying",
)
@click.option(
    "-D",
    "--save-diff",
    type=str,
    metavar="<filename>",
    default=None,
    help="Save diffs to a file instead of updating in place",
)
@click.option(
    "-u",
    "--update-deps",
    is_flag=True,
    default=False,
    help="Automatically update all linked modules and subworkflows without asking for confirmation",
)
def command_subworkflows_update(
    ctx,
    subworkflow,
    directory,
    force,
    prompt,
    sha,
    install_all,
    preview,
    save_diff,
    update_deps,
    limit_output,
):
    """
    Update DSL2 subworkflow within a pipeline.
    """
    subworkflows_update(
        ctx, subworkflow, directory, force, prompt, sha, install_all, preview, save_diff, update_deps, limit_output
    )


# nf-core test-dataset subcommands
@nf_core_cli.group(aliases=["tds"])
@click.pass_context
def test_datasets(ctx):
    """
    Commands to manage nf-core test datasets.
    """
    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)


# nf-core test-dataset search
@test_datasets.command("search", short_help="Search files in the nf-core/test-datasets repository")
@click.pass_context
@click.option("-b", "--branch", type=str, help="Branch in the test-datasets repository to reduce search to")
@click.option(
    "-p",
    "--generate-nf-path",
    is_flag=True,
    default=False,
    help="Auto-generate a file path for use in nextflow code based on the branch and query result",
)
@click.option(
    "-u",
    "--generate-dl-url",
    is_flag=True,
    default=False,
    help="Auto-generate a github url for downloading the test data file based on the branch and query result. Only applicable when not using -p / --generate-nf-path",
)
@click.argument("query", required=False)
def command_test_dataset_search(ctx, branch, generate_nf_path, generate_dl_url, query):
    """
    Search files filtered by QUERY on a specified branch in the nf-core/test-datasets repository.
    If no QUERY is given or QUERY is ambiguous, an auto-completion form is shown.
    """
    test_datasets_search(ctx, branch, generate_nf_path, generate_dl_url, query)


# nf-core test-dataset search
@test_datasets.command("list")
@click.pass_context
@click.option("-b", "--branch", type=str, help="Branch in the test-datasets repository to reduce search to")
@click.option(
    "-p",
    "--generate-nf-path",
    is_flag=True,
    default=False,
    help="Auto-generate a file path for use in nextflow code based on the branch and query result",
)
@click.option(
    "-u",
    "--generate-dl-url",
    is_flag=True,
    default=False,
    help="Auto-generate a github url for downloading the test data file based on the branch and query result. Only applicable when not using -p / --generate-nf-path",
)
def command_test_dataset_list_remote(ctx, branch, generate_nf_path, generate_dl_url):
    """
    List files on a specified branch in the nf-core/test-datasets repository.
    """
    test_datasets_list_remote(ctx, branch, generate_nf_path, generate_dl_url)


# nf-core test-datasets list-branches
@test_datasets.command("list-branches")
@click.pass_context
def command_test_datasets_list_branches(ctx):
    """
    List remote branches with test data in the nf-core/test-dataset repository.
    """
    test_datasets_list_branches(ctx)


## DEPRECATED commands since v3.0.0


# nf-core schema subcommands (deprecated)
@nf_core_cli.group(deprecated=True, hidden=True)
def schema():
    """
    Use `nf-core pipelines schema <command>` instead.
    """
    pass


# nf-core schema validate (deprecated)
@schema.command("validate", deprecated=True)
@click.argument("pipeline", required=True, metavar="<pipeline name>")
@click.argument("params", type=click.Path(exists=True), required=True, metavar="<JSON params file>")
def command_schema_validate(pipeline, params):
    """
    Use `nf-core pipelines schema validate` instead.
    """
    log.warning(
        "The `[magenta]nf-core schema validate[/]` command is deprecated. Use `[magenta]nf-core pipelines schema validate[/]` instead."
    )
    pipelines_schema_validate(pipeline, params)


# nf-core schema build (deprecated)
@schema.command("build", deprecated=True)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "--no-prompts",
    is_flag=True,
    help="Do not confirm changes, just update parameters and exit",
)
@click.option(
    "--web-only",
    is_flag=True,
    help="Skip building using Nextflow config, just launch the web tool",
)
@click.option(
    "--url",
    type=str,
    default="https://oldsite.nf-co.re/pipeline_schema_builder",
    help="Customise the builder URL (for development work)",
)
def command_schema_build(directory, no_prompts, web_only, url):
    """
    Use `nf-core pipelines schema build` instead.
    """
    log.warning(
        "The `[magenta]nf-core schema build[/]` command is deprecated. Use `[magenta]nf-core pipelines schema build[/]` instead."
    )
    pipelines_schema_build(directory, no_prompts, web_only, url)


# nf-core schema lint (deprecated)
@schema.command("lint", deprecated=True)
@click.argument(
    "schema_path",
    type=click.Path(exists=True),
    default="nextflow_schema.json",
    metavar="<pipeline schema>",
)
def command_schema_lint(schema_path):
    """
    Use `nf-core pipelines schema lint` instead.
    """
    log.warning(
        "The `[magenta]nf-core schema lint[/]` command is deprecated. Use `[magenta]nf-core pipelines schema lint[/]` instead."
    )
    pipelines_schema_lint(schema_path)


# nf-core schema docs (deprecated)
@schema.command("docs", deprecated=True)
@click.argument(
    "schema_path",
    type=click.Path(exists=True),
    default="nextflow_schema.json",
    required=False,
    metavar="<pipeline schema>",
)
@click.option(
    "-o",
    "--output",
    type=str,
    metavar="<filename>",
    help="Output filename. Defaults to standard out.",
)
@click.option(
    "-x",
    "--format",
    type=click.Choice(["markdown", "html"]),
    default="markdown",
    help="Format to output docs in.",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-c",
    "--columns",
    type=str,
    metavar="<columns_list>",
    help="CSV list of columns to include in the parameter tables (parameter,description,type,default,required,hidden)",
    default="parameter,description,type,default,required,hidden",
)
def command_schema_docs(schema_path, output, format, force, columns):
    """
    Use `nf-core pipelines schema docs` instead.
    """
    log.warning(
        "The `[magenta]nf-core schema docs[/]` command is deprecated. Use `[magenta]nf-core pipelines schema docs[/]` instead."
    )
    pipelines_schema_docs(schema_path, output, format, force, columns)


# nf-core create-logo (deprecated)
@nf_core_cli.command("create-logo", deprecated=True, hidden=True)
@click.argument("logo-text", metavar="<logo_text>")
@click.option("-d", "--dir", "directory", type=click.Path(), default=".", help="Directory to save the logo in.")
@click.option(
    "-n",
    "--name",
    type=str,
    help="Name of the output file (with or without '.png' suffix).",
)
@click.option(
    "--theme",
    type=click.Choice(["light", "dark"]),
    default="light",
    help="Theme for the logo.",
    show_default=True,
)
@click.option(
    "--width",
    type=int,
    default=2300,
    help="Width of the logo in pixels.",
    show_default=True,
)
@click.option(
    "--format",
    type=click.Choice(["png", "svg"]),
    default="png",
    help="Image format of the logo, either PNG or SVG.",
    show_default=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite any files if they already exist",
)
def command_create_logo(logo_text, directory, name, theme, width, format, force):
    """
    Use `nf-core pipelines create-logo` instead.
    """
    log.warning(
        "The `[magenta]nf-core create-logo[/]` command is deprecated. Use `[magenta]nf-core pipeliness create-logo[/]` instead."
    )
    pipelines_create_logo(logo_text, directory, name, theme, width, format, force)


# nf-core sync (deprecated)
@nf_core_cli.command("sync", hidden=True, deprecated=True)
@click.pass_context
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-b",
    "--from-branch",
    type=str,
    help="The git branch to use to fetch workflow variables.",
)
@click.option(
    "-p",
    "--pull-request",
    is_flag=True,
    default=False,
    help="Make a GitHub pull-request with the changes.",
)
@click.option(
    "--force_pr",
    is_flag=True,
    default=False,
    help="Force the creation of a pull-request, even if there are no changes.",
)
@click.option("-g", "--github-repository", type=str, help="GitHub PR: target repository.")
@click.option("-u", "--username", type=str, help="GitHub PR: auth username.")
@click.option("-t", "--template-yaml", help="Pass a YAML file to customize the template")
def command_sync(ctx, directory, from_branch, pull_request, github_repository, username, template_yaml, force_pr):
    """
    Use `nf-core pipelines sync` instead.
    """
    log.warning(
        "The `[magenta]nf-core sync[/]` command is deprecated. Use `[magenta]nf-core pipelines sync[/]` instead."
    )
    pipelines_sync(ctx, directory, from_branch, pull_request, github_repository, username, template_yaml, force_pr)


# nf-core bump-version (deprecated)
@nf_core_cli.command("bump-version", hidden=True, deprecated=True)
@click.pass_context
@click.argument("new_version", default="")
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory. [dim]\[default: current working directory][/]",
)
@click.option(
    "-n",
    "--nextflow",
    is_flag=True,
    default=False,
    help="Bump required nextflow version instead of pipeline version",
)
def command_bump_version(ctx, new_version, directory, nextflow):
    """
    Use `nf-core pipelines bump-version` instead.
    """
    log.warning(
        "The `[magenta]nf-core bump-version[/]` command is deprecated. Use `[magenta]nf-core pipelines bump-version[/]` instead."
    )
    pipelines_bump_version(ctx, new_version, directory, nextflow)


# nf-core list (deprecated)
@nf_core_cli.command("list", deprecated=True, hidden=True)
@click.argument("keywords", required=False, nargs=-1, metavar="<filter keywords>")
@click.option(
    "-s",
    "--sort",
    type=click.Choice(["release", "pulled", "name", "stars"]),
    default="release",
    help="How to sort listed pipelines",
)
@click.option("--json", is_flag=True, default=False, help="Print full output as JSON")
@click.option("--show-archived", is_flag=True, default=False, help="Print archived workflows")
@click.pass_context
def command_list(ctx, keywords, sort, json, show_archived):
    """
    Use `nf-core pipelines list` instead.
    """
    log.warning(
        "The `[magenta]nf-core list[/]` command is deprecated. Use `[magenta]nf-core pipelines list[/]` instead."
    )
    pipelines_list(ctx, keywords, sort, json, show_archived)


# nf-core launch (deprecated)
@nf_core_cli.command("launch", deprecated=True, hidden=True)
@click.argument("pipeline", required=False, metavar="<pipeline name>")
@click.option("-r", "--revision", help="Release/branch/SHA of the project to run (if remote)")
@click.option("-i", "--id", help="ID for web-gui launch parameter set")
@click.option(
    "-c",
    "--command-only",
    is_flag=True,
    default=False,
    help="Create Nextflow command with params (no params file)",
)
@click.option(
    "-o",
    "--params-out",
    type=click.Path(),
    default=os.path.join(os.getcwd(), "nf-params.json"),
    help="Path to save run parameters file",
)
@click.option(
    "-p",
    "--params-in",
    type=click.Path(exists=True),
    help="Set of input run params to use from a previous run",
)
@click.option(
    "-a",
    "--save-all",
    is_flag=True,
    default=False,
    help="Save all parameters, even if unchanged from default",
)
@click.option(
    "-x",
    "--show-hidden",
    is_flag=True,
    default=False,
    help="Show hidden params which don't normally need changing",
)
@click.option(
    "-u",
    "--url",
    type=str,
    default="https://nf-co.re/launch",
    help="Customise the builder URL (for development work)",
)
@click.pass_context
def command_launch(
    ctx,
    pipeline,
    id,
    revision,
    command_only,
    params_in,
    params_out,
    save_all,
    show_hidden,
    url,
):
    """
    Use `nf-core pipelines launch` instead.
    """
    log.warning(
        "The `[magenta]nf-core launch[/]` command is deprecated. Use `[magenta]nf-core pipelines launch[/]` instead."
    )
    pipelines_launch(ctx, pipeline, id, revision, command_only, params_in, params_out, save_all, show_hidden, url)


# nf-core create-params-file (deprecated)
@nf_core_cli.command("create-params-file", deprecated=True, hidden=True)
@click.argument("pipeline", required=False, metavar="<pipeline name>")
@click.option("-r", "--revision", help="Release/branch/SHA of the pipeline (if remote)")
@click.option(
    "-o",
    "--output",
    type=str,
    default="nf-params.yml",
    metavar="<filename>",
    help="Output filename. Defaults to `nf-params.yml`.",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-x",
    "--show-hidden",
    is_flag=True,
    default=False,
    help="Show hidden params which don't normally need changing",
)
def command_create_params_file(pipeline, revision, output, force, show_hidden):
    """
    Use `nf-core pipelines create-params-file` instead.
    """
    log.warning(
        "The `[magenta]nf-core create-params-file[/]` command is deprecated. Use `[magenta]nf-core pipelines create-params-file[/]` instead."
    )
    pipelines_create_params_file(pipeline, revision, output, force, show_hidden)


# nf-core download (deprecated)
@nf_core_cli.command("download", deprecated=True, hidden=True)
@click.argument("pipeline", required=False, metavar="<pipeline name>")
@click.option(
    "-r",
    "--revision",
    multiple=True,
    help="Pipeline release to download. Multiple invocations are possible, e.g. `-r 1.1 -r 1.2`",
)
@click.option("-o", "--outdir", type=str, help="Output directory")
@click.option(
    "-x",
    "--compress",
    type=click.Choice(["tar.gz", "tar.bz2", "zip", "none"]),
    help="Archive compression type",
)
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite existing files")
@click.option(
    "-t",
    "--tower",
    is_flag=True,
    default=False,
    hidden=True,
    help="Download for Seqera Platform. DEPRECATED: Please use `--platform` instead.",
)
@click.option(
    "--platform",
    is_flag=True,
    default=False,
    help="Download for Seqera Platform (formerly Nextflow Tower)",
)
@click.option(
    "-d",
    "--download-configuration",
    is_flag=True,
    default=False,
    help="Include configuration profiles in download. Not available with `--platform`",
)
@click.option(
    "--tag",
    multiple=True,
    help="Add custom alias tags to `--platform` downloads. For example, `--tag \"3.10=validated\"` adds the custom 'validated' tag to the 3.10 release.",
)
@click.option(
    "-s",
    "--container-system",
    type=click.Choice(["none", "singularity", "docker"]),
    help="Download container images of required software.",
)
@click.option(
    "-l",
    "--container-library",
    multiple=True,
    help="Container registry/library or mirror to pull images from.",
)
@click.option(
    "-u",
    "--container-cache-utilisation",
    type=click.Choice(["amend", "copy", "remote"]),
    help="Utilise a `singularity.cacheDir` in the download process, if applicable.",
)
@click.option(
    "-i",
    "--container-cache-index",
    type=str,
    help="List of images already available in a remote `singularity.cacheDir`.",
)
@click.option(
    "-p",
    "--parallel-downloads",
    type=int,
    default=4,
    help="Number of parallel image downloads",
)
@click.pass_context
def command_download(
    ctx,
    pipeline,
    revision,
    outdir,
    compress,
    force,
    tower,
    platform,
    download_configuration,
    tag,
    container_system,
    container_library,
    container_cache_utilisation,
    container_cache_index,
    parallel_downloads,
):
    """
    Use `nf-core pipelines download` instead.
    """
    log.warning(
        "The `[magenta]nf-core download[/]` command is deprecated. Use `[magenta]nf-core pipelines download[/]` instead."
    )
    pipelines_download(
        ctx,
        pipeline,
        revision,
        outdir,
        compress,
        force,
        platform or tower,
        download_configuration,
        tag,
        container_system,
        container_library,
        container_cache_utilisation,
        container_cache_index,
        parallel_downloads,
    )


# nf-core lint (deprecated)
@nf_core_cli.command("lint", hidden=True, deprecated=True)
@click.option(
    "-d",
    "--dir",
    "directory",
    type=click.Path(exists=True),
    default=".",
    help=r"Pipeline directory [dim]\[default: current working directory][/]",
)
@click.option(
    "--release",
    is_flag=True,
    default=Path(os.environ.get("GITHUB_REF", "").strip(" '\"")).parent.name in ["master", "main"]
    and os.environ.get("GITHUB_REPOSITORY", "").startswith("nf-core/")
    and not os.environ.get("GITHUB_REPOSITORY", "") == "nf-core/tools",
    help="Execute additional checks for release-ready workflows.",
)
@click.option(
    "-f",
    "--fix",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Attempt to automatically fix specified lint test",
)
@click.option(
    "-k",
    "--key",
    type=str,
    metavar="<test>",
    multiple=True,
    help="Run only these lint tests",
)
@click.option("-p", "--show-passed", is_flag=True, help="Show passing tests on the command line")
@click.option("-i", "--fail-ignored", is_flag=True, help="Convert ignored tests to failures")
@click.option("-w", "--fail-warned", is_flag=True, help="Convert warn tests to failures")
@click.option(
    "--markdown",
    type=str,
    metavar="<filename>",
    help="File to write linting results to (Markdown)",
)
@click.option(
    "--json",
    type=str,
    metavar="<filename>",
    help="File to write linting results to (JSON)",
)
@click.option(
    "--sort-by",
    type=click.Choice(["module", "test"]),
    default="test",
    help="Sort lint output by module or test name.",
    show_default=True,
)
@click.pass_context
def command_lint(
    ctx,
    directory,
    release,
    fix,
    key,
    show_passed,
    fail_ignored,
    fail_warned,
    markdown,
    json,
    sort_by,
):
    """
    Use `nf-core pipelines lint` instead.
    """
    log.warning(
        "The `[magenta]nf-core lint[/]` command is deprecated. Use `[magenta]nf-core pipelines lint[/]` instead."
    )
    pipelines_lint(ctx, directory, release, fix, key, show_passed, fail_ignored, fail_warned, markdown, json, sort_by)


# nf-core create (deprecated)
@nf_core_cli.command("create", hidden=True, deprecated=True)
@click.option(
    "-n",
    "--name",
    type=str,
    help="The name of your new pipeline",
)
@click.option("-d", "--description", type=str, help="A short description of your pipeline")
@click.option("-a", "--author", type=str, help="Name of the main author(s)")
@click.option("--version", type=str, default="1.0.0dev", help="The initial version number to use")
@click.option("-f", "--force", is_flag=True, default=False, help="Overwrite output directory if it already exists")
@click.option("-o", "--outdir", help="Output directory for new pipeline (default: pipeline name)")
@click.option("-t", "--template-yaml", help="Pass a YAML file to customize the template")
@click.option("--plain", is_flag=True, help="Use the standard nf-core template")
@click.option(
    "--organisation",
    type=str,
    default="nf-core",
    help="The name of the GitHub organisation where the pipeline will be hosted (default: nf-core)",
)
@click.pass_context
def command_create(ctx, name, description, author, version, force, outdir, template_yaml, plain, organisation):
    """
    Use `nf-core pipelines create` instead.
    """
    log.warning(
        "The `[magenta]nf-core create[/]` command is deprecated. Use `[magenta]nf-core pipelines create[/]` instead."
    )
    pipelines_create(ctx, name, description, author, version, force, outdir, template_yaml, organisation)


# Main script is being run - launch the CLI
if __name__ == "__main__":
    run_nf_core()
