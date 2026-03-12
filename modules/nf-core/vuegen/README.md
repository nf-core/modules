# VueGen Module

## Dockerimage

The dockerimage is needed as both quarto and the TinyTeX distribution are installed
inside the container to make vuegen work properly. At runtime especially apptainer
does not support installing packages on the fly to read-only directories.

```bash
#!/bin/bash
set -euo pipefail

ORG="quay.io/nf-core"
VERSION="v0.5.1"

# 1. Build & push the base image (nextflow version)
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  --push \
  -t ${ORG}/vuegen:${VERSION} \
  -f Dockerfile .
```
