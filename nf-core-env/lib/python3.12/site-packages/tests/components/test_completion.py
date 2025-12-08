import pytest

from nf_core.components.components_completion import autocomplete_modules, autocomplete_subworkflows

from ..utils import GITLAB_NFTEST_BRANCH, GITLAB_URL


class DummyParam:
    # Minimal mock object for Click parameter (not used in the function)
    pass


class DummyCtx:
    def __init__(self, obj=None, params=None):
        self.obj = obj
        self.params = params if params is not None else {}


def test_autocomplete_modules():
    ctx = DummyCtx(
        obj={
            "modules_repo_url": GITLAB_URL,
            "modules_repo_branch": GITLAB_NFTEST_BRANCH,
            "modules_repo_no_pull": True,
        }
    )
    param = DummyParam()
    completions = autocomplete_modules(ctx, param, "samt")

    values = [c.value for c in completions]
    assert "samtools/stats" in values
    assert "samtools/idxstats" in values
    assert "fastqc" not in values


def test_autocomplete_modules_missing_argument(capfd):
    ctx = DummyCtx()
    param = DummyParam()

    with pytest.raises(TypeError) as exc_info:
        autocomplete_modules(ctx, param)  # Missing 'incomplete' argument

    assert "missing 1 required positional argument" in str(exc_info.value)


def test_autocomplete_subworkflows():
    ctx = DummyCtx(
        obj={
            "modules_repo_url": GITLAB_URL,
            "modules_repo_branch": GITLAB_NFTEST_BRANCH,
            "modules_repo_no_pull": True,
        }
    )
    param = DummyParam()
    completions = autocomplete_subworkflows(ctx, param, "bam_stats")

    values = [c.value for c in completions]
    assert "bam_stats_samtools" in values
    assert "bam_sort_stats_samtools" not in values


def test_autocomplete_subworkflows_missing_argument():
    ctx = DummyCtx()
    param = DummyParam()

    with pytest.raises(TypeError) as exc_info:
        autocomplete_subworkflows(ctx, param)  # Missing 'incomplete' argument

    assert "missing 1 required positional argument" in str(exc_info.value)
