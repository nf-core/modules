import logging
import subprocess
from pathlib import Path

import yaml
from rich.console import Console
from rich.logging import RichHandler

# Configure rich console
console = Console()

# Configure logging with rich
logging.basicConfig(
    level="DEBUG",  # Set to DEBUG to capture all messages
    format="%(message)s",
    datefmt="[%X]",
    handlers=[
        RichHandler(console=console, rich_tracebacks=True, show_time=False, level="INFO"),
    ],
)
logger = logging.getLogger("rich")


def get_wave_image(env_file):
    logger.debug(f"Requesting Wave image for {env_file}")
    try:
        result = subprocess.run(
            ["wave", "--conda-file", env_file, "--await"], capture_output=True, text=True, check=True
        )
        logger.debug(f"Successfully obtained Wave image for {env_file}")
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running Wave CLI for {env_file}: {e}")
        return None


def update_meta_yml(meta_file, docker_image):
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
        containers_index = next(i for i, line in enumerate(lines) if line.strip() == "containers:")
        indent = len(lines[containers_index + 1]) - len(lines[containers_index + 1].lstrip())
        lines.insert(containers_index + 1, " " * indent + f"docker_x86: {docker_image}")
        content = "\n".join(lines) + "\n"

    # Write the updated content back to the file
    with open(meta_file, "w") as f:
        f.write(content)

    logger.debug(f"Successfully updated {meta_file}")


def process_module_directory(directory, modname):
    env_file = directory / "environment.yml"
    meta_file = directory / "meta.yml"

    if not env_file.exists():
        logger.debug(f"Skipping {directory.name}: environment.yml not found")
        return
    if not meta_file.exists():
        logger.debug(f"Skipping {directory.name}: meta.yml not found")
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


def process_modules():
    modules_dir = Path("modules/nf-core")
    logger.debug(f"Processing modules in {modules_dir}")

    for module_dir in modules_dir.iterdir():
        if not module_dir.is_dir():
            continue

        # Check if the module directory itself contains environment.yml and meta.yml
        if (module_dir / "environment.yml").exists() and (module_dir / "meta.yml").exists():
            process_module_directory(module_dir, module_dir.name)
        else:
            # If not, check for subdirectories
            for subdir in module_dir.iterdir():
                if subdir.is_dir():
                    logger.debug(f"Processing subdirectory: {subdir.name}")
                    process_module_directory(subdir, f"{module_dir.name}/{subdir.name}")


if __name__ == "__main__":
    console.rule("[bold blue]Starting module processing[/bold blue]")
    process_modules()
    console.rule("[bold blue]Finished processing all modules[/bold blue]")
