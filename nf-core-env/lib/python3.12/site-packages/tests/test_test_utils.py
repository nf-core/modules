from pathlib import Path

from .utils import with_temporary_file, with_temporary_folder


def test_with_temporary_file():
    @with_temporary_file
    def tmp_file_exists(tmp_file):
        assert Path(tmp_file.name).exists()

    tmp_file_exists()


def test_does_not_exist_after():
    tmp_file = with_temporary_file(lambda x: x.name)()
    assert not Path(tmp_file).exists()


def test_with_temporary_folder():
    @with_temporary_folder
    def tmp_folder_exists(tmp_folder):
        assert Path(tmp_folder).exists()

    tmp_folder_exists()


def test_tmp_folder_does_not_exist_after():
    tmp_folder = with_temporary_folder(lambda x: x)()
    assert not Path(tmp_folder).exists()
