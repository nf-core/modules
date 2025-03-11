# Dockerfile to create container with Cell Ranger Arc v2.0.2 and bcl2fastq v2.20.0
# Push to quay.io/nf-core/cellranger-arc-mkfastq:<VER>

FROM continuumio/miniconda3:4.8.2
LABEL authors="Regina Reynolds, Florian Heyl, Gisela Gabernet <gisela.gabernet@gmail.com>" \
    description="Docker image containing bcl2fastq2 and Cell Ranger Arc"
# Disclaimer: this container is not provided nor supported by Illumina or 10x Genomics.

# Install procps and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    cpio \
    procps \
    rpm2cpio \
    unzip \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded bcl2fastq2 and cellranger-arc file
ENV BCL2FASTQ2_VER=v2-20-0-linux-x86-64 \
    ARC_VER=2.0.2
COPY bcl2fastq2-$BCL2FASTQ2_VER.zip /tmp/bcl2fastq2-$BCL2FASTQ2_VER.zip
COPY cellranger-arc-$ARC_VER.tar.gz /opt/cellranger-arc-$ARC_VER.tar.gz

# Install bcl2fastq2
RUN \
    cd /tmp && \
    unzip bcl2fastq2-$BCL2FASTQ2_VER.zip && \
    mv *.rpm bcl2fastq2-$BCL2FASTQ2_VER.rpm && \
    rpm2cpio ./bcl2fastq2-$BCL2FASTQ2_VER.rpm | cpio -idmv && \
    export PATH=/tmp/usr/local/bin/:$PATH && \
    mv /tmp/usr/local/bin/bcl2fastq /usr/bin/bcl2fastq && \
    rm -rf bcl2fastq2-$BCL2FASTQ2_VER.*

# Install cellranger-arc
RUN \
    cd /opt && \
    tar -xzvf cellranger-arc-$ARC_VER.tar.gz && \
    export PATH=/opt/cellranger-arc-$ARC_VER:$PATH && \
    ln -s /opt/cellranger-arc-$ARC_VER/cellranger-arc /usr/bin/cellranger-arc && \
    rm -rf /opt/cellranger-arc-$ARC_VER.tar.gz
