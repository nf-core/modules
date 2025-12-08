"""Context manager for temporarily changing folder."""

import os

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


class FolderContext(object):
    """Context manager for temporarily changing directory."""

    def __init__(self, folder):
        """
        Store the previous working path to restore upon exit.

        :param str folder: Path to set as new working directory
        """
        if not os.path.isdir(folder):
            raise ValueError("Requested temp entry to non-folder: {}".format(folder))
        self._prevdir = os.getcwd()
        self._currdir = folder

    def __enter__(self):
        """Make the working directory switch."""
        os.chdir(self._currdir)

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Switch back to the previous working directory."""
        if not os.path.isdir(self._prevdir):
            raise RuntimeError(
                "Return path is no longer a directory: {}".format(self._prevdir)
            )
        os.chdir(self._prevdir)
