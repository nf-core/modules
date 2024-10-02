#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "httpx",
# ]
# ///

import logging

logger = logging.getLogger(__name__)

image_url = "oras://community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:aa20de1f1b5ddb30"

wave_api_url = "https://wave.seqera.io"
url = f"{wave_api_url}/v1alpha1/inspect"
logger.info(f"calling image inspect at {url} for image url {image_url}")
data = {"containerimage": image_url}

# if platform_pat:
#     data["toweraccesstoken"] = platform_pat
# else:
# TODO
logger.warning("'platform_pat' not set, no auth to wave back end")

with httpx.client() as client:
    response = client.post(
        url=url,
        json=data,
        headers={"content-type": "application/json"},
    )

    if response.is_success:
        data = response.json()
        layers = data.get("container", {}).get("manifest", {}).get("layers", [])
        is_singularity = len(layers) == 1 and layers[0].get("mediatype", "").endswith(".sif")
        if not is_singularity:
            raise httpexception(status_code=400, detail="not a singularity image")
        if "digest" not in layers[0]:
            raise httpexception(status_code=400, detail="no 'digest' in first layer found")
        digest = layers[0]["digest"].replace("sha256:", "")
        container_url = (
            f"https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/{digest[:2]}/{digest}/data"
        )
        print(container_url)

    else:
        raise httpexception(status_code=response.status_code, detail=response.text)
