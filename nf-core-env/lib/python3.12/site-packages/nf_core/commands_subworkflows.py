import logging
import sys

import rich

from nf_core.utils import rich_force_colors

log = logging.getLogger(__name__)

stdout = rich.console.Console(force_terminal=rich_force_colors())


def subworkflows_create(ctx, subworkflow, directory, author, force, migrate_pytest):
    """
    Create a new subworkflow from the nf-core template.

    If the specified directory is a pipeline, this function creates a file called
    'subworkflows/local/<subworkflow_name>.nf'

    If the specified directory is a clone of nf-core/modules, it creates or modifies files
    in 'subworkflows/' and 'tests/subworkflows'
    """
    from nf_core.subworkflows import SubworkflowCreate

    # Run function
    try:
        subworkflow_create = SubworkflowCreate(directory, subworkflow, author, force, migrate_pytest)
        subworkflow_create.create()
    except UserWarning as e:
        log.critical(e)
        sys.exit(1)
    except LookupError as e:
        log.error(e)
        sys.exit(1)


def subworkflows_test(ctx, subworkflow, directory, no_prompts, update, once, profile, migrate_pytest):
    """
    Run nf-test for a subworkflow.

    Given the name of a subworkflow, runs the nf-test command to test the subworkflow and generate snapshots.
    """
    from nf_core.components.components_test import ComponentsTest

    if migrate_pytest:
        subworkflows_create(ctx, subworkflow, directory, None, False, True)
    try:
        sw_tester = ComponentsTest(
            component_type="subworkflows",
            component_name=subworkflow,
            directory=directory,
            no_prompts=no_prompts,
            update=update,
            once=once,
            remote_url=ctx.obj["modules_repo_url"],
            branch=ctx.obj["modules_repo_branch"],
            verbose=ctx.obj["verbose"],
            profile=profile,
        )
        sw_tester.run()
    except (UserWarning, LookupError) as e:
        log.critical(e)
        sys.exit(1)


def subworkflows_list_remote(ctx, keywords, json):
    """
    List subworkflows in a remote GitHub repo [dim i](e.g [link=https://github.com/nf-core/modules]nf-core/modules[/])[/].
    """
    from nf_core.subworkflows import SubworkflowList

    try:
        subworkflow_list = SubworkflowList(
            ".",
            True,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )

        stdout.print(subworkflow_list.list_components(keywords, json))
    except (UserWarning, LookupError) as e:
        log.critical(e)
        sys.exit(1)


def subworkflows_list_local(ctx, keywords, json, directory):  # pylint: disable=redefined-builtin
    """
    List subworkflows installed locally in a pipeline
    """
    from nf_core.subworkflows import SubworkflowList

    try:
        subworkflow_list = SubworkflowList(
            directory,
            False,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )
        stdout.print(subworkflow_list.list_components(keywords, json))
    except (UserWarning, LookupError) as e:
        log.error(e)
        sys.exit(1)


def subworkflows_lint(ctx, subworkflow, directory, registry, key, all, fail_warned, local, passed, sort_by, fix):
    """
    Lint one or more subworkflows in a directory.

    Checks DSL2 subworkflow code against nf-core guidelines to ensure
    that all subworkflows follow the same standards.

    Test subworkflows within a pipeline or a clone of the
    nf-core/modules repository.
    """
    from nf_core.components.lint import LintExceptionError
    from nf_core.subworkflows import SubworkflowLint

    try:
        subworkflow_lint = SubworkflowLint(
            directory,
            fail_warned=fail_warned,
            fix=fix,
            registry=ctx.params["registry"],
            remote_url=ctx.obj["modules_repo_url"],
            branch=ctx.obj["modules_repo_branch"],
            no_pull=ctx.obj["modules_repo_no_pull"],
            hide_progress=ctx.obj["hide_progress"],
        )
        subworkflow_lint.lint(
            subworkflow=subworkflow,
            registry=registry,
            key=key,
            all_subworkflows=all,
            print_results=True,
            local=local,
            show_passed=passed,
            sort_by=sort_by,
        )
        if len(subworkflow_lint.failed) > 0:
            sys.exit(1)
    except LintExceptionError as e:
        log.error(e)
        sys.exit(1)
    except (UserWarning, LookupError) as e:
        log.critical(e)
        sys.exit(1)


def subworkflows_info(ctx, subworkflow, directory):
    """
    Show developer usage information about a given subworkflow.

    Parses information from a subworkflow's [i]meta.yml[/] and renders help
    on the command line. A handy equivalent to searching the
    [link=https://nf-co.re/modules]nf-core website[/].

    If run from a pipeline and a local copy of the subworkflow is found, the command
    will print this usage info.
    If not, usage from the remote subworkflows repo will be shown.
    """
    from nf_core.subworkflows import SubworkflowInfo

    try:
        subworkflow_info = SubworkflowInfo(
            directory,
            subworkflow,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )
        stdout.print(subworkflow_info.get_component_info())
    except (UserWarning, LookupError) as e:
        log.error(e)
        sys.exit(1)


def subworkflows_install(ctx, subworkflow, directory, prompt, force, sha):
    """
    Install DSL2 subworkflow within a pipeline.

    Fetches and installs subworkflow files from a remote repo e.g. nf-core/modules.
    """
    from nf_core.subworkflows import SubworkflowInstall

    try:
        subworkflow_install = SubworkflowInstall(
            directory,
            force,
            prompt,
            sha,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )
        exit_status = subworkflow_install.install(subworkflow)
        if not exit_status:
            sys.exit(1)
    except (UserWarning, LookupError) as e:
        log.error(e)
        sys.exit(1)


def subworkflows_remove(ctx, directory, subworkflow):
    """
    Remove a subworkflow from a pipeline.
    """
    from nf_core.subworkflows import SubworkflowRemove

    try:
        module_remove = SubworkflowRemove(
            directory,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
        )
        module_remove.remove(subworkflow)
    except (UserWarning, LookupError) as e:
        log.critical(e)
        sys.exit(1)


def subworkflows_update(
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

    Fetches and updates subworkflow files from a remote repo e.g. nf-core/modules.
    """
    from nf_core.subworkflows import SubworkflowUpdate

    try:
        subworkflow_install = SubworkflowUpdate(
            directory,
            force,
            prompt,
            sha,
            install_all,
            preview,
            save_diff,
            update_deps,
            ctx.obj["modules_repo_url"],
            ctx.obj["modules_repo_branch"],
            ctx.obj["modules_repo_no_pull"],
            limit_output,
        )
        exit_status = subworkflow_install.update(subworkflow)
        if not exit_status and install_all:
            sys.exit(1)
    except (UserWarning, LookupError) as e:
        log.error(e)
        sys.exit(1)
