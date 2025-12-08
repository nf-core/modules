import shutil

import git
import pytest

import nf_core.pipelines.lint_utils

JSON_WITH_SYNTAX_ERROR = "{'a':1, 1}"
JSON_MALFORMED = "{'a':1}"
JSON_FORMATTED = '{ "a": 1 }\n'

WHICH_PRE_COMMIT = shutil.which("pre-commit")


@pytest.fixture()
def temp_git_repo(tmp_path_factory):
    tmp_git_dir = tmp_path_factory.mktemp("tmp_git_dir")
    repo = git.Repo.init(tmp_git_dir)
    return tmp_git_dir, repo


@pytest.fixture(name="formatted_json")
def git_dir_with_json(temp_git_repo):
    tmp_git_dir, repo = temp_git_repo
    file = tmp_git_dir / "formatted.json"
    with open(file, "w", encoding="utf-8") as f:
        f.write(JSON_FORMATTED)
    repo.git.add(file)
    return file


@pytest.fixture(name="malformed_json")
def git_dir_with_json_malformed(temp_git_repo):
    tmp_git_dir, repo = temp_git_repo
    file = tmp_git_dir / "malformed.json"
    with open(file, "w", encoding="utf-8") as f:
        f.write(JSON_MALFORMED)
    repo.git.add(file)
    return file


@pytest.fixture(name="syntax_error_json")
def git_dir_with_json_syntax_error(temp_git_repo):
    tmp_git_dir, repo = temp_git_repo
    file = tmp_git_dir / "synthax-error.json"
    with open(file, "w", encoding="utf-8") as f:
        f.write(JSON_WITH_SYNTAX_ERROR)
    repo.git.add(file)
    return file


def test_run_prettier_on_formatted_file(formatted_json):
    nf_core.pipelines.lint_utils.run_prettier_on_file(formatted_json)
    assert formatted_json.read_text() == JSON_FORMATTED


def test_run_prettier_on_malformed_file(malformed_json):
    nf_core.pipelines.lint_utils.run_prettier_on_file(malformed_json)
    assert malformed_json.read_text() == JSON_FORMATTED


def test_run_prettier_on_syntax_error_file(syntax_error_json, caplog):
    nf_core.pipelines.lint_utils.run_prettier_on_file(syntax_error_json)
    expected_critical_log = "SyntaxError: Unexpected token (1:10)"
    assert expected_critical_log in caplog.text
