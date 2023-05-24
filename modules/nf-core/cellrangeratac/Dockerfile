# Dockerfile to create container with Cell Ranger ATAC v2.1.0
# Push to nfcore/cellrangeratac:<VER>

FROM continuumio/miniconda3:4.8.2
LABEL authors="Gisela Gabernet <gisela.gabernet@gmail.com>, Florian Heyl" \
    description="Docker image containing Cell Ranger"
# Disclaimer: this container is not provided nor supported by Illumina or 10x Genomics.

# Install procps and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    cpio \
    procps \
    rpm2cpio \
    unzip \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded cellrangeratac file
ENV CELLRANGERATAC_VER=2.1.0
COPY cellranger-atac-$CELLRANGERATAC_VER.tar.gz /opt/cellranger-atac-$CELLRANGERATAC_VER.tar.gz

# Install cellrangeratac
RUN \
    cd /opt && \
    tar -xzvf cellranger-atac-$CELLRANGERATAC_VER.tar.gz && \
    export PATH=/opt/cellranger-atac-$CELLRANGERATAC_VER:$PATH && \
    ln -s /opt/cellranger-atac-$CELLRANGERATAC_VER/cellranger-atac /usr/bin/cellranger-atac && \
    rm -rf /opt/cellranger-atac-$CELLRANGERATAC_VER.tar.gz
