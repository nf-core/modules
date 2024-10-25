#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "httpx",
# ]
# ///

import logging

import httpx

logger = logging.getLogger(__name__)

image_url = "oras://community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:aa20de1f1b5ddb30"

if image_url.startswith("oras://"):
    image_url = image_url.replace("oras://", "")

wave_api_url = "https://wave.seqera.io"
url = f"{wave_api_url}/v1alpha1/inspect"

# if platform_pat:
#     data["toweraccesstoken"] = platform_pat
# else:
# TODO
logger.warning("'platform_pat' not set, no auth to wave back end")

try:
    logger.info(f"calling image inspect at {url} for image url {image_url}")
    response = httpx.post(
        url=url,
        json={"containerImage": image_url},
        headers={"content-type": "application/json"},
    )

    data = response.json()
    logger.debug(data)
    layers = data.get("container", {}).get("manifest", {}).get("layers", [])
    is_singularity = len(layers) == 1 and layers[0].get("mediaType", "").endswith(".sif")
    if not is_singularity:
        print(layers)
        raise ValueError("not a singularity image")
    if "digest" not in layers[0]:
        print(layers)
        raise ValueError("no 'digest' in first layer found")

    digest = layers[0]["digest"].replace("sha256:", "")
    container_url = f"https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/{digest[:2]}/{digest}/data"
    print(container_url)

except httpx.RequestError as exc:
    print(f"An error occurred while requesting {exc.request.url!r}.")
    print("No singularity image for you")
