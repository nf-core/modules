# Dockerfile to create a container with dorado
# Push to quay.io/nf-core/dorado:<VER>
#
# The vendor image (docker.io/nanoporetech/dorado) is 9.4 GB compressed because
# it also ships a second copy of the CUDA runtime and 4.3 GB of basecalling
# models. This image carries only the official release tarball, which bundles the
# libraries dorado links against, so it serves both the CPU and the GPU modules.
# Models are downloaded at run time or supplied as a module input.

FROM ubuntu:22.04 AS builder

ARG DORADO_VER=1.4.0
ARG DORADO_SHA256=ea15250d911c6b6efdb21cc142008b2be8a5c781c8b0a4593b82f27b07385731

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL -o /dorado.tar.gz https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VER}-linux-x64.tar.gz \
    && echo "${DORADO_SHA256}  /dorado.tar.gz" | sha256sum -c - \
    && mkdir -p /opt/dorado \
    && tar -xzf /dorado.tar.gz -C /opt/dorado --strip-components 1 \
    && rm /dorado.tar.gz

# multistage so that the tarball does not stay behind in a layer
FROM ubuntu:22.04

LABEL authors="Samuel Ahuno" \
    description="Docker image containing dorado"
# Disclaimer: this container is not provided nor supported by Oxford Nanopore Technologies.

# 'ps' is needed by Nextflow to collect task metrics; dorado needs CA
# certificates, and falls back to curl, to download models
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/dorado /opt/dorado

ENV PATH="/opt/dorado/bin:${PATH}"
