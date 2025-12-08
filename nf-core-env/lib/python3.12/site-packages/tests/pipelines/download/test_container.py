"""Tests for the download subcommand of nf-core tools"""

import unittest

import pytest
import rich.progress_bar
import rich.table
import rich.text

from nf_core.pipelines.download.container_fetcher import ContainerProgress


class ContainerProgressTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    #
    # Test for 'utils..add/update_main_task'
    #
    def test_download_progress_main_task(self):
        with ContainerProgress() as progress:
            # No task initially
            assert progress.tasks == []

            # Add a task, it should be there
            task_id = progress.add_main_task(total=42)
            assert task_id == 0
            assert len(progress.tasks) == 1
            assert progress.task_ids[0] == task_id
            assert progress.tasks[0].total == 42

            # Add another task, there should now be two
            other_task_id = progress.add_task("Another task", total=28)
            assert other_task_id == 1
            assert len(progress.tasks) == 2
            assert progress.task_ids[1] == other_task_id
            assert progress.tasks[1].total == 28

            progress.update_main_task(total=35)
            assert progress.tasks[0].total == 35
            assert progress.tasks[1].total == 28

    #
    # Test for 'utils.DownloadProgress.sub_task'
    #
    def test_download_progress_sub_task(self):
        with ContainerProgress() as progress:
            # No task initially
            assert progress.tasks == []

            # Add a sub-task, it should be there
            with progress.sub_task("Sub-task", total=42) as sub_task_id:
                assert sub_task_id == 0
                assert len(progress.tasks) == 1
                assert progress.task_ids[0] == sub_task_id
                assert progress.tasks[0].total == 42

            # The sub-task should be gone now
            assert progress.tasks == []

            # Add another sub-task, this time that raises an exception
            with pytest.raises(ValueError):
                with progress.sub_task("Sub-task", total=28) as sub_task_id:
                    assert sub_task_id == 1
                    assert len(progress.tasks) == 1
                    assert progress.task_ids[0] == sub_task_id
                    assert progress.tasks[0].total == 28
                    raise ValueError("This is a test error")

            # The sub-task should also be gone now
            assert progress.tasks == []

    #
    # Test for 'utils.DownloadProgress.get_renderables'
    #
    def test_download_progress_renderables(self):
        # Test the "summary" progress type
        with ContainerProgress() as progress:
            assert progress.tasks == []
            progress.add_task("Task 1", progress_type="summary", total=42, completed=11)
            assert len(progress.tasks) == 1

            renderable = progress.get_renderable()
            assert isinstance(renderable, rich.console.Group), type(renderable)

            assert len(renderable.renderables) == 1
            table = renderable.renderables[0]
            assert isinstance(table, rich.table.Table)

            assert isinstance(table.columns[0]._cells[0], str)
            assert table.columns[0]._cells[0] == "[magenta]Task 1"

            assert isinstance(table.columns[1]._cells[0], rich.progress_bar.ProgressBar)
            assert table.columns[1]._cells[0].completed == 11
            assert table.columns[1]._cells[0].total == 42

            assert isinstance(table.columns[2]._cells[0], str)
            assert table.columns[2]._cells[0] == "[progress.percentage] 26%"

            assert isinstance(table.columns[3]._cells[0], str)
            assert table.columns[3]._cells[0] == "•"

            assert isinstance(table.columns[4]._cells[0], str)
            assert table.columns[4]._cells[0] == "[green]11/42 tasks completed"


class ContainerTest(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def use_caplog(self, caplog):
        self._caplog = caplog

    @property
    def logged_levels(self) -> list[str]:
        return [record.levelname for record in self._caplog.records]

    @property
    def logged_messages(self) -> list[str]:
        return [record.message for record in self._caplog.records]

    def __contains__(self, item: str) -> bool:
        """Allows to check for log messages easily using the in operator inside a test:
        assert 'my log message' in self
        """
        return any(record.message == item for record in self._caplog.records if self._caplog)

    pass
