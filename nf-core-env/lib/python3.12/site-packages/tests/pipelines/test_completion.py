import pytest

from nf_core.pipelines.list import autocomplete_pipelines


class DummyParam:
    # Minimal mock object for Click parameter
    pass


class DummyCtx:
    def __init__(self, obj=None, params=None):
        self.obj = obj
        self.params = params if params is not None else {}


def test_autocomplete_pipelines():
    ctx = DummyCtx()
    param = DummyParam()
    completions = autocomplete_pipelines(ctx, param, "nf-core/bac")

    values = [c.value for c in completions]
    print(values)  # For debugging purposes

    assert "nf-core/bacass" in values
    assert "nf-core/bactmap" in values
    assert "nf-core/abotyper" not in values


def test_autocomplete_pipelines_missing_argument(capfd):
    ctx = DummyCtx()
    param = DummyParam()

    with pytest.raises(TypeError) as exc_info:
        autocomplete_pipelines(ctx, param)  # Missing 'incomplete' argument

    assert "missing 1 required positional argument" in str(exc_info.value)
