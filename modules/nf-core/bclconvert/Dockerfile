# syntax=docker/dockerfile:1
# Dockerfile to create container with bcl-convert
# Push to nfcore/bclconvert:<VER>

FROM debian:bullseye-slim
ARG BCLCONVERT_VERSION="4.3.6"
LABEL org.opencontainers.image.authors="Matthias De Smet <matthias.desmet@ugent.be>"
LABEL org.opencontainers.image.description="Docker image containing bcl-convert"
LABEL org.opencontainers.image.version="$BCLCONVERT_VERSION"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/LICENSE"

# Disclaimer: this container is not provided nor supported by Illumina
# 'ps' command is need by some nextflow executions to collect system stats
# Install procps and clean apt cache
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    rpm2cpio \
    cpio \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname
# Install bcl-convert
ADD bcl-convert-4.3.6-2.el8.x86_64.rpm bcl-convert.rpm
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN rpm2cpio bcl-convert.rpm | cpio -idmv \
    && rm bcl-convert.rpm
RUN rm -rf /var/log/bcl-convert && ln -sfT /tmp /var/log/bcl-convert
