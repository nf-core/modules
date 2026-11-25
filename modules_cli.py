import argparse
import hashlib
import logging
import os
import subprocess
from pathlib import Path

import yaml
from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn


console = Console()


def configure_logging():
    """Configure rich logging output."""
    logging.basicConfig(
        level="DEBUG",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[
            RichHandler(
                console=console,
                rich_tracebacks=True,
                show_time=False,
                level="INFO",
            ),
        ],
    )
    return logging.getLogger("rich")


logger = configure_logging()


# ---------------------------------------------------------------------------
# Conda environment collection (from get_conda_envs.py)
# ---------------------------------------------------------------------------

def find_and_normalize_yaml(module: str | None = None):
    """
    Find all environment.yml files under modules/nf-core (optionally limited to
    a single module) and write de-duplicated, normalised copies to conda_envs/.

    If `module` is provided, it should be a path relative to modules/nf-core,
    e.g. "tcoffee/align".
    """
    base_dir = Path("modules/nf-core")
    if module:
        base_dir = base_dir / module

    output_dir = Path("conda_envs")
    output_dir.mkdir(exist_ok=True)

    processed_hashes: set[str] = set()
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

        for root, _, files in os.walk(base_dir):
            for file in files:
                if file == "environment.yml":
                    total_files += 1
                    file_path = Path(root) / file
                    was_written = normalize_and_save_yaml(
                        file_path, output_dir, processed_hashes
                    )
                    if was_written:
                        written_files += 1
                    else:
                        skipped_files += 1

                    progress.update(task, advance=1)

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


def normalize_and_save_yaml(
    file_path: Path, output_dir: Path, processed_hashes: set[str]
) -> bool:
    with open(file_path) as f:
        data = yaml.safe_load(f)

    normalized_yaml = yaml.dump(data, sort_keys=False)
    yaml_hash = hashlib.md5(normalized_yaml.encode()).hexdigest()

    if yaml_hash not in processed_hashes:
        output_file = output_dir / f"{yaml_hash}.yml"
        with open(output_file, "w") as f:
            f.write(normalized_yaml)
        processed_hashes.add(yaml_hash)
        console.print(
            f"[green]✅   Saved:[/] {file_path} [dim](location: {output_file})[/]"
        )
        return True

    console.print(f"[yellow]⏭️ Skipped:[/] {file_path}")
    return False


# ---------------------------------------------------------------------------
# Wave image generation and meta.yml update (from get_wave_images.py)
# ---------------------------------------------------------------------------

def get_wave_image(env_file: str) -> str | None:
    logger.debug(f"Requesting Wave image for {env_file}")
    try:
        # Always use --freeze when generating images
        result = subprocess.run(
            ["wave", "--conda-file", env_file, "--await", "--freeze"],
            capture_output=True,
            text=True,
            check=True,
        )
        logger.debug(f"Successfully obtained Wave image for {env_file}")
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running Wave CLI for {env_file}: {e}")
        return None


def update_meta_yml(meta_file: Path, docker_image: str) -> None:
    logger.debug(f"Updating {meta_file} with Docker image")

    # Read the existing content
    with open(meta_file) as f:
        content = f.read()
        meta_data = yaml.safe_load(content) or {}

    # Check if 'containers' key exists
    if "containers" not in meta_data:
        # If not, add it with proper indentation
        content += f"\ncontainers:\n  docker_x86: {docker_image}\n"
    else:
        # If it exists, add docker_x86 under it
        lines = content.splitlines()
        containers_index = next(
            i for i, line in enumerate(lines) if line.strip() == "containers:"
        )
        indent = len(lines[containers_index + 1]) - len(
            lines[containers_index + 1].lstrip()
        )
        lines.insert(containers_index + 1, " " * indent + f"docker_x86: {docker_image}")
        content = "\n".join(lines) + "\n"

    # Write the updated content back to the file
    with open(meta_file, "w") as f:
        f.write(content)

    logger.debug(f"Successfully updated {meta_file}")


def process_module_directory(directory: Path, modname: str) -> None:
    env_file = directory / "environment.yml"
    meta_file = directory / "meta.yml"

    if not env_file.exists():
        logger.debug(f"Skipping {directory}: environment.yml not found")
        return
    if not meta_file.exists():
        logger.debug(f"Skipping {directory}: meta.yml not found")
        return

    logger.debug(f"Checking {meta_file}")
    with open(meta_file) as f:
        meta_data = yaml.safe_load(f) or {}

    if "containers" in meta_data and "docker_x86" in meta_data["containers"]:
        logger.debug(f"Skipped {meta_file} (already contains docker_x86)")
        return

    logger.debug(f"docker_x86 not found in {meta_file}, proceeding with Wave CLI")
    logger.info(f"Processing: {modname}")
    docker_image = get_wave_image(str(env_file))
    if docker_image:
        update_meta_yml(meta_file, docker_image)
        logger.debug(f"Got image: {docker_image}")
    else:
        logger.error(f"Failed to update {meta_file}: No Docker image obtained")


def process_modules(module: str | None = None) -> None:
    """
    Walk modules/nf-core and populate docker_x86 using Wave images.

    If `module` is provided, only that module (e.g. "tcoffee/align") is processed.
    """
    modules_dir = Path("modules/nf-core")
    logger.debug(f"Processing modules in {modules_dir}")

    # If a specific module path is requested, restrict traversal to that path
    if module:
        target = modules_dir / module
        if not target.exists():
            logger.error(f"Requested module '{module}' not found under {modules_dir}")
            return

        # If the target itself is a module directory with the files, process it
        if (target / "environment.yml").exists() and (target / "meta.yml").exists():
            process_module_directory(target, module)
            return

        # Otherwise, search its subdirectories
        for subdir in target.iterdir():
            if subdir.is_dir():
                modname = f"{module}/{subdir.name}"
                logger.debug(f"Processing subdirectory: {modname}")
                process_module_directory(subdir, modname)
        return

    # No specific module requested: process all modules
    for module_dir in modules_dir.iterdir():
        if not module_dir.is_dir():
            continue

        # Check if the module directory itself contains environment.yml and meta.yml
        if (module_dir / "environment.yml").exists() and (
            module_dir / "meta.yml"
        ).exists():
            process_module_directory(module_dir, module_dir.name)
        else:
            # If not, check for subdirectories
            for subdir in module_dir.iterdir():
                if subdir.is_dir():
                    modname = f"{module_dir.name}/{subdir.name}"
                    logger.debug(f"Processing subdirectory: {modname}")
                    process_module_directory(subdir, modname)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect unique Conda environments and generate Wave images to "
            "populate docker_x86 entries in nf-core module meta.yml files."
        )
    )
    parser.add_argument(
        "-m",
        "--module",
        help=(
            "Only process a single nf-core module, "
            "e.g. 'tcoffee/align'. By default all modules are processed."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()

    console.rule("[bold magenta]Starting Conda environment processing[/bold magenta]")
    find_and_normalize_yaml(module=args.module)

    console.rule("[bold blue]Starting Wave image generation[/bold blue]")
    process_modules(module=args.module)
    console.rule("[bold blue]Finished processing modules[/bold blue]")


if __name__ == "__main__":
    main()


