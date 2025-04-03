# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.10-jammy
ARG SPACERANGER_VER="3.1.3"
# NOTE you need to go fill out the form to get a new download link
ARG SPACERANGER_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.3.tar.gz?Expires=1739251752&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=ohLEz5jR~HNzWrUz4pd35VGlJP~g7ywmVDaD29BoCdUvio46HWEfhIuUT19D~2JoFlYsbaLQ2d~qIuBB-d9WmX9bdmoQ1vzc9ozikwufXhuJhvrvxIXAZ4kUreoGWX7b51emM7nq8kJxTL9CmTKQZXRe~DuQjP6q6z0XxDUxeLhdPRg3ACfj~jFgMlNfwiJwBRHDH3yM5Hrtz~gNmxZ9O0fl~KNwlzZh7Mtm9GG-GyXYOtM5KAPQjv69JpM2XiZoZFpaVbPq5uqGfr~B3GC1dY25MGqVgXavc~uQFVcsZiyF5BewbeFzYTKu0SiW4WdO21K5L6CdkJ1M5AflXZmpsw__" \"
ARG SPACERANGER_SHA256="cee7e94942cf72af76e28ead0024992642f81f10b89006bf93756c1e08528522"
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
ADD --checksum=sha256:$SPACERANGER_SHA256 \
    $SPACERANGER_URL \
    /opt/spaceranger-$SPACERANGER_VER.tar.gz

# Install spaceranger
RUN tar -xzf /opt/spaceranger-$SPACERANGER_VER.tar.gz -C /opt/ \
    && rm /opt/spaceranger-$SPACERANGER_VER.tar.gz \
    && ln -s /opt/spaceranger-$SPACERANGER_VER/spaceranger /usr/bin/spaceranger
ENV PATH="/opt/spaceranger-$SPACERANGER_VER:$PATH"
