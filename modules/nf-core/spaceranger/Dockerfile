# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.10-jammy
ARG SPACERANGER_VER="4.0.1"
ARG SPACERANGER_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-4.0.1.tar.gz?Expires=1756263067&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ZuXyhILJpqvX67FsKT4HuAW4JAW0pkLCDphqWCpqe553ujNcB9r13EmZU-yS7Wwmrfr7iPPYZsBGs2CEHQMjgJa0EhnB3kAIgTu33mtnBlZscmjbPx~AkgmvCJxrF06hXfU7nPIG0nQH2TbwAV2XAY67~-IHn8xrLZhEYdIFedIcsbpXJ0RMsb~JDyCaqqykIRnDsS85fAJ-nO7zaiv53iUsN4PgvVwF9vBJslLQzH-pA44ihpA0C4ibHxXIztdKbo72Gu-udibDvmiXBoM3Sb7E1-W6VCp8z0R9Yzwpv1POjYel7bVVc8bZd1Nj3vNNamr2T0fyP1XHnHcbDCZlEg__" \"
ARG SPACERANGER_SHA256="eacb4a3e7a342424a57ac757fce95d536c1e0099fd8fef8d249f8afb105d6d91"
LABEL org.opencontainers.image.authors="Gisela Gabernet <gisela.gabernet@gmail.com>,Edmund Miller <edmund@nf-co.re>"
LABEL org.opencontainers.image.description="Docker image containing spaceranger"
LABEL org.opencontainers.image.version="$SPACERANGER_VER"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/spaceranger/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://www.10xgenomics.com/support/software/space-ranger/downloads/eula"
# Disclaimer: this container is not provided nor supported by 10x Genomics.

# Install procps and curl and clean apt cache
RUN micromamba install -y -n base conda-forge::procps-ng \
    && micromamba clean -a -y
USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

# NOTE Docker doesn't support md5sum and that's what spaceranger gives us --checksum=md5sum:949088b1e2cc5bd03079b60a65b41d39
# https://docs.docker.com/reference/dockerfile/#add---checksum
# hadolint ignore=DL3020
ADD --checksum=sha256:$SPACERANGER_SHA256 \
    $SPACERANGER_URL \
    /opt/spaceranger-$SPACERANGER_VER.tar.gz

# Install spaceranger
RUN tar -xzf /opt/spaceranger-$SPACERANGER_VER.tar.gz -C /opt/ \
    && rm /opt/spaceranger-$SPACERANGER_VER.tar.gz \
    && ln -s /opt/spaceranger-$SPACERANGER_VER/spaceranger /usr/bin/spaceranger
ENV PATH="/opt/spaceranger-$SPACERANGER_VER:$PATH"
