# syntax=docker/dockerfile:1
# Dockerfile to create container with bcl-convert
# Push to nfcore/bclconvert:<VER>

# Build stage: unpack the RPM file
FROM debian:bullseye-slim AS build
ARG BCLCONVERT_VERSION="4.4.6"

# Install tools needed to extract RPM
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    rpm2cpio \
    cpio \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Extract bcl-convert from RPM
COPY bcl-convert-4.4.6-2.el8.x86_64.rpm bcl-convert.rpm
RUN rpm2cpio bcl-convert.rpm | cpio -idmv \
    && rm bcl-convert.rpm

# Final stage: minimal image with only necessary components
FROM debian:bullseye-slim
ARG BCLCONVERT_VERSION="4.4.6"
LABEL org.opencontainers.image.description="Docker image containing bcl-convert"
LABEL org.opencontainers.image.version="$BCLCONVERT_VERSION"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://github.com/nf-core/modules/blob/master/modules/nf-core/bclconvert/LICENSE"

# Disclaimer: this container is not provided nor supported by Illumina
# Install procps (for 'ps' command needed by Nextflow) and rclone (for parallel file transfers)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    rclone \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname

# Copy bcl-convert executable from build stage
COPY --from=build /usr/bin/bcl-convert /usr/bin/bcl-convert
