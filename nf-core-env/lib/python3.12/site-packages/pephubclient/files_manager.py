import os
from contextlib import suppress
from pathlib import Path

import pandas
import yaml
import zipfile

from pephubclient.exceptions import PEPExistsError


class FilesManager:
    @staticmethod
    def save_jwt_data_to_file(path: str, jwt_data: str) -> None:
        """
        Save jwt to provided path
        """
        Path(os.path.dirname(path)).mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write(jwt_data)

    @staticmethod
    def load_jwt_data_from_file(path: str) -> str:
        """
        Open the file with username and ID and load this data.
        """
        with suppress(FileNotFoundError):
            with open(path, "r") as f:
                return f.read()

    @staticmethod
    def create_project_folder(
        parent_path: str,
        folder_name: str,
    ) -> str:
        """
        Create new project folder

        :param parent_path: parent path to create folder in
        :param folder_name: folder name
        :return: folder_path
        """
        if parent_path:
            if not Path(parent_path).exists():
                raise OSError(
                    f"Parent path does not exist. Provided path: {parent_path}"
                )
        folder_path = os.path.join(parent_path or os.getcwd(), folder_name)
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        return folder_path

    @staticmethod
    def save_yaml(config: dict, full_path: str, not_force: bool = False):
        FilesManager.check_writable(path=full_path, force=not not_force)
        with open(full_path, "w") as outfile:
            yaml.dump(config, outfile, default_flow_style=False)

    @staticmethod
    def save_pandas(df: pandas.DataFrame, full_path: str, not_force: bool = False):
        FilesManager.check_writable(path=full_path, force=not not_force)
        df.to_csv(full_path, index=False)

    @staticmethod
    def file_exists(full_path: str) -> bool:
        return os.path.isfile(full_path)

    @staticmethod
    def delete_file_if_exists(filename: str) -> None:
        with suppress(FileNotFoundError):
            os.remove(filename)
            print(
                f"\033[38;5;11m{f'File was deleted successfully -> {filename}'}\033[0m"
            )

    @staticmethod
    def check_writable(path: str, force: bool = True):
        if not force and os.path.isfile(path):
            raise PEPExistsError(f"File already exists and won't be updated: {path}")

    @staticmethod
    def save_zip_file(files_dict: dict, file_path: str, force: bool = False) -> None:
        """
        Save zip file with provided files as dict.

        :param files_dict: dict with files to save. e.g. {"file1.txt": "file1 content"}
        :param file_path: filename to save zip file to
        :param force: overwrite file if exists
        :return: None
        """
        FilesManager.check_writable(path=file_path, force=force)
        with zipfile.ZipFile(
            file_path, mode="w", compression=zipfile.ZIP_DEFLATED
        ) as zf:
            for name, res in files_dict.items():
                zf.writestr(name, str.encode(res))
