FROM mambaorg/micromamba:1.4-focal

LABEL authors="Nicolas Vannieuwkerke <nicolas.vannieuwkerke@ugent.be>" \
    description="Docker image containing AnnotSV"

ENV VERSION=3.3.2
ENV ANNOTSV=/tmp/AnnotSV-${VERSION}
ENV PATH="${PATH}:${ANNOTSV}/bin"

# Install mamba dependencies
RUN micromamba install -y --name base -c conda-forge -c bioconda -c defaults \
    git \
    make \
    wget \
    coreutils \
    tk=8.6.12 \
    bedtools=2.30 \
    bcftools=1.16 \
    openjdk=8 \
    unzip \
    python=3.8 \
    natsort=7.1.1 \
    pandas=1.5.3 \
    pyfaidx=0.7.2.1 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v$VERSION.tar.gz && \
    tar -xvf v$VERSION.tar.gz -C . && \
    cd AnnotSV* && \
    make PREFIX=${ANNOTSV} install
