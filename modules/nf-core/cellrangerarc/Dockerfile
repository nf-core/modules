# Dockerfile to create container with Cell Ranger v2.0.2
# Push to quay.io/nf-core/cellranger-arc:<VER>

FROM continuumio/miniconda3:4.8.2
LABEL authors="Gisela Gabernet <gisela.gabernet@gmail.com>, Florian Heyl" \
    description="Docker image containing Cell Ranger Arc"
# Disclaimer: this container is not provided nor supported by Illumina or 10x Genomics.

# Install procps and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    cpio \
    procps \
    rpm2cpio \
    unzip \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded cellranger-arc file
ENV CELLRANGER_ARC_VER=2.0.2
COPY cellranger-arc-$CELLRANGER_ARC_VER.tar.gz /opt/cellranger-arc-$CELLRANGER_ARC_VER.tar.gz

# Install cellranger-arc
RUN \
    cd /opt && \
    tar -xzvf cellranger-arc-$CELLRANGER_ARC_VER.tar.gz && \
    export PATH=/opt/cellranger-arc-$CELLRANGER_ARC_VER:$PATH && \
    ln -s /opt/cellranger-arc-$CELLRANGER_ARC_VER/cellranger-arc /usr/bin/cellranger-arc && \
    rm -rf /opt/cellranger-arc-$CELLRANGER_ARC_VER.tar.gz
