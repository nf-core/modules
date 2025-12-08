import os
import subprocess
import unittest
from pathlib import Path

import pytest

from nf_core.pipelines.download.utils import (
    DownloadError,
    intermediate_file,
)

from ...utils import with_temporary_folder


class DownloadUtilsTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    #
    # Test for 'utils.intermediate_file'
    #
    @with_temporary_folder
    def test_intermediate_file(self, outdir):
        outdir = Path(outdir)
        # Code that doesn't fail. The file shall exist

        # Directly write to the file, as in download_image
        output_path = outdir / "testfile1"
        with intermediate_file(output_path) as tmp:
            tmp_path = Path(tmp.name)
            tmp.write(b"Hello, World!")

        assert output_path.exists()
        assert os.path.getsize(output_path) == 13
        assert not tmp_path.exists()

        # Run an external command as in pull_image
        output_path = outdir / "testfile2"
        with intermediate_file(output_path) as tmp:
            tmp_path = Path(tmp.name)
            subprocess.check_call([f"echo 'Hello, World!' > {tmp_path}"], shell=True)

        assert (output_path).exists()
        assert os.path.getsize(output_path) == 14  # Extra \n !
        assert not (tmp_path).exists()

        # Code that fails. The file shall not exist

        # Directly write to the file and raise an exception
        output_path = outdir / "testfile3"
        with pytest.raises(ValueError):
            with intermediate_file(output_path) as tmp:
                tmp_path = Path(tmp.name)
                tmp.write(b"Hello, World!")
                raise ValueError("This is a test error")

        assert not (output_path).exists()
        assert not (tmp_path).exists()

        # Run an external command and raise an exception
        output_path = outdir / "testfile4"
        with pytest.raises(subprocess.CalledProcessError):
            with intermediate_file(output_path) as tmp:
                tmp_path = Path(tmp.name)
                subprocess.check_call([f"echo 'Hello, World!' > {tmp_path}"], shell=True)
                subprocess.check_call(["ls", "/dummy"])

        assert not (output_path).exists()
        assert not (tmp_path).exists()

        # Test for invalid output paths
        with pytest.raises(DownloadError):
            with intermediate_file(outdir) as tmp:
                pass

        output_path = outdir / "testfile5"
        os.symlink("/dummy", output_path)
        with pytest.raises(DownloadError):
            with intermediate_file(output_path) as tmp:
                pass
