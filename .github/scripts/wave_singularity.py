#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "requests",
#     "rich",
#     "rich_click"
# ]
# ///

import logging

import requests
import rich_click as click
from rich.logging import RichHandler

# Replace the basic logger setup with rich logging
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        RichHandler(
            rich_tracebacks=True,
            show_time=False,
            markup=True,
        )
    ],
)
logger = logging.getLogger(__name__)
click.rich_click.SHOW_ARGUMENTS = True


@click.command()
@click.option(
    "--platform-pat",
    envvar="SEQERA_ACCESS_TOKEN",
    show_envvar=True,
    help="Platform authentication token",
)
@click.argument("image_name")
def main(image_name, platform_pat):
    """Script to return a HTTPS Singularity image URL from Seqera Containers."""

    if image_name.startswith("oras://"):
        image_name = image_name.replace("oras://", "")

    wave_api_url_base = "https://wave.seqera.io"
    wave_api_url = f"{wave_api_url_base}/v1alpha1/inspect"
    container_url_base = "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256"

    request_data = {"containerImage": image_name}
    if platform_pat:
        request_data["toweraccesstoken"] = platform_pat
    else:
        logger.debug("'--platform-pat' / '$TOWER_ACCESS_TOKEN' not set, no auth to wave back end")

    try:
        logger.debug(f"Calling image inspect at {wave_api_url} for image url '{image_name}'")
        response = requests.post(
            url=wave_api_url,
            json=request_data,
            headers={"content-type": "application/json"},
        )

        data = response.json()
        logger.debug(data)
        layers = data.get("container", {}).get("manifest", {}).get("layers", [])
        is_singularity = len(layers) == 1 and layers[0].get("mediaType", "").endswith(".sif")
        logger.debug(layers)
        if not is_singularity:
            raise ValueError("Not a singularity image")
        if "digest" not in layers[0]:
            raise ValueError("no 'digest' in first layer found")

        digest = layers[0]["digest"].replace("sha256:", "")
        container_url = f"{container_url_base}/{digest[:2]}/{digest}/data"
        print(container_url)

    except requests.RequestException as exc:
        raise ValueError(f"An error occurred while requesting {wave_api_url}\n {exc}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        logger.error(f"[red]{exc}[/red]")
        exit(1)
