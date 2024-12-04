# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.10-noble

ARG BCL2FASTQ_VERSION="2.20.0"
ARG BCL2FASTQ_URL="https://support.illumina.com/softwaredownload.html?assetId=82660c4c-f46c-4743-8566-2437755e4329&assetDetails=bcl2fastq2-v2-20-0-tar.zip"
ARG BCL2FASTQ_SHA256="2566b24f29829b39f3add112a674990b1c54ae2fbe7ccb50a4c7dce9ccf152e6"

LABEL org.opencontainers.image.authors="Matthias De Smet <matthias.desmet@ugent.be>, Edmund Miller <edmund@nf-co.re>"
LABEL org.opencontainers.image.description="Docker image containing bcl2fastq"
LABEL org.opencontainers.image.version="$BCL2FASTQ_VERSION"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/bcl2fastq/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf"
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
# Install bcl2fastq
ADD bcl2fastq-2.20.0.422-1.el8.x86_64.rpm bcl2fastq.rpm
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN rpm2cpio bcl2fastq.rpm | cpio -idmv \
    && rm bcl2fastq.rpm
RUN rm -rf /var/log/bcl2fastq && ln -sfT /tmp /var/log/bcl2fastq


# Install procps and curl and clean apt cache
RUN micromamba install -y -n base conda-forge::procps-ng \
    && micromamba clean -a -y
USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

# NOTE Docker doesn't support md5sum and that's what spaceranger gives us --checksum=md5sum:949088b1e2cc5bd03079b60a65b41d39
# https://docs.docker.com/reference/dockerfile/#add---checksum
ADD --checksum=sha256:$SPACERANGER_SHA256 \
    $SPACERANGER_URL \
    /opt/spaceranger-$SPACERANGER_VER.tar.gz

# Install spaceranger
RUN tar -xzf /opt/spaceranger-$SPACERANGER_VER.tar.gz -C /opt/ \
    && rm /opt/spaceranger-$SPACERANGER_VER.tar.gz \
    && ln -s /opt/spaceranger-$SPACERANGER_VER/spaceranger /usr/bin/spaceranger
ENV PATH="/opt/spaceranger-$SPACERANGER_VER:$PATH"
