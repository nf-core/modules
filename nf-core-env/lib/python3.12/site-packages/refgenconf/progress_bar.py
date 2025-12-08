from datetime import timedelta

from rich.progress import ProgressColumn, filesize
from rich.text import Text


class _DownloadColumn(ProgressColumn):
    """Renders file size downloaded and total, e.g. '0.5/2.3 GB'."""

    @staticmethod
    def render(task):
        """Calculate common unit for completed and total."""
        completed = int(task.completed)
        total = int(task.total)
        unit, suffix = filesize.pick_unit_and_suffix(
            total, ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"], 1024
        )
        completed_ratio = completed / unit
        total_ratio = total / unit
        precision = 0 if unit == 1 else 1
        completed_str = f"{completed_ratio:,.{precision}f}"
        total_str = f"{total_ratio:,.{precision}f}"
        download_status = f"{completed_str}/{total_str} {suffix}"
        download_text = Text(download_status, style="[bright_white]")
        return download_text


class _TransferSpeedColumn(ProgressColumn):
    """Renders human readable transfer speed."""

    @staticmethod
    def render(task):
        """Show data transfer speed."""
        speed = task.speed
        if speed is None:
            return Text("?", style="[bright_white]")
        data_speed = filesize.decimal(int(speed))
        return Text(f"{data_speed}/s", style="[bright_white]")


class _TimeRemainingColumn(ProgressColumn):
    """Renders estimated time remaining."""

    # Only refresh twice a second to prevent jitter
    max_refresh = 1

    @staticmethod
    def render(task):
        """Show time remaining."""
        remaining = task.time_remaining
        if remaining is None:
            return Text("-:--:--", style="[bright_white]")
        remaining_delta = timedelta(seconds=int(remaining))
        return Text(str(remaining_delta), style="[bright_white]")
