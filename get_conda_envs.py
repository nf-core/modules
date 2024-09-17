import hashlib
import os
from pathlib import Path

import yaml
from rich.console import Console
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn

console = Console()


def find_and_normalize_yaml():
    base_dir = Path("modules/nf-core")
    output_dir = Path("conda_envs")
    output_dir.mkdir(exist_ok=True)

    processed_hashes = set()
    total_files = 0
    written_files = 0
    skipped_files = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        console=console,
    ) as progress:
        task = progress.add_task("[green]Processing files...", total=None)

        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if file == "environment.yml":
                    total_files += 1
                    file_path = Path(root) / file
                    was_written = normalize_and_save_yaml(file_path, output_dir, processed_hashes, progress, task)
                    if was_written:
                        written_files += 1
                    else:
                        skipped_files += 1

    console.print(
        Panel.fit(
            f"[bold green]✨ Summary ✨[/]\n\n"
            f"📁 Total environment.yml files found: [cyan]{total_files}[/]\n"
            f"📝 Unique files written: [cyan]{written_files}[/]\n"
            f"🔄 Duplicate files skipped: [cyan]{skipped_files}[/]",
            title="Conda Environment Processing Results",
            border_style="bold blue",
        )
    )


def normalize_and_save_yaml(file_path, output_dir, processed_hashes, progress, task):
    with open(file_path) as f:
        data = yaml.safe_load(f)

    normalized_yaml = yaml.dump(data, sort_keys=False)
    yaml_hash = hashlib.md5(normalized_yaml.encode()).hexdigest()

    if yaml_hash not in processed_hashes:
        output_file = output_dir / f"{yaml_hash}.yml"
        with open(output_file, "w") as f:
            f.write(normalized_yaml)
        processed_hashes.add(yaml_hash)
        progress.console.print(f"[green]✅   Saved:[/] {file_path} [dim](location: {output_file})[/]")
        return True
    else:
        progress.console.print(f"[yellow]⏭️ Skipped:[/] {file_path}")
        return False

    progress.update(task, advance=1)


if __name__ == "__main__":
    console.print("[bold magenta]Starting Conda Environment Processing...[/]\n")
    find_and_normalize_yaml()
