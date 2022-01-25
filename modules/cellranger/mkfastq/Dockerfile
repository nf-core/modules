# Dockerfile to create container with Cell Ranger v6.1.2 and bcl2fastq v2.20.0
# Push to nfcore/cellrangermkfastq:<VER>

FROM continuumio/miniconda3:4.8.2
LABEL authors="Regina Reynolds, Gisela Gabernet <gisela.gabernet@gmail.com>" \
    description="Docker image containing bcl2fastq2 and Cell Ranger"
# Disclaimer: this container is not provided nor supported by Illumina or 10x Genomics.

# Install procps and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    cpio \
    procps \
    rpm2cpio \
    unzip \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded bcl2fastq2 and cellranger file
ENV BCL2FASTQ2_VER=v2-20-0-linux-x86-64 \
    CELLRANGER_VER=6.1.2
COPY bcl2fastq2-$BCL2FASTQ2_VER.zip /tmp/bcl2fastq2-$BCL2FASTQ2_VER.zip
COPY cellranger-$CELLRANGER_VER.tar.gz /opt/cellranger-$CELLRANGER_VER.tar.gz

# Install bcl2fastq2
RUN \
    cd /tmp && \
    unzip bcl2fastq2-$BCL2FASTQ2_VER.zip && \
    mv *.rpm bcl2fastq2-$BCL2FASTQ2_VER.rpm && \
    rpm2cpio ./bcl2fastq2-$BCL2FASTQ2_VER.rpm | cpio -idmv && \
    export PATH=/tmp/usr/local/bin/:$PATH && \
    ln -s /tmp/usr/local/bin/bcl2fastq /usr/bin/bcl2fastq && \
    rm -rf bcl2fastq2-$BCL2FASTQ2_VER.*

# Install cellranger
RUN \
    cd /opt && \
    tar -xzvf cellranger-$CELLRANGER_VER.tar.gz && \
    export PATH=/opt/cellranger-$CELLRANGER_VER:$PATH && \
    ln -s /opt/cellranger-$CELLRANGER_VER/cellranger /usr/bin/cellranger && \
    rm -rf /opt/cellranger-$CELLRANGER_VER.tar.gz
