# syntax=docker/dockerfile:1

FROM mambaorg/micromamba:1.5.10-noble
ARG SPACERANGER_VER="3.1.2"
# NOTE you need to go fill out the form to get a new download link
ARG SPACERANGER_URL="https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.2.tar.gz?Expires=1732608367&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=S4jcGCs5H0wLFfREhSc~AfuPIOXE~PW96kX9o2EyxZnmn1goUecgfEWTj67fL1DnZJBIa06kFdUdNpsMn2ustiSWOfXEgjlYQnxIcNnQIiEBGoQTHWphvx3LEQ6wtZnkWS80P6IcE0HJkIsgy04t6Sohih5cxY4jgytYsrAfZDYr5G3KKFwTfCKmhzMaXqW635yPbyQ8xEcQHK0QwviAx8-EFq-PE8UzC4QgUKi2MW-ivcfZkSDSfF8C3s7SgwDXIGIWv52mzeszenxMjN4KrWQotZ7ZpktzI0Vfpz0dNC17dQeDQUHj4LuNYbdh3RqsPKtqu3wjCe2Q7KiyoWnmaw__" \"
ARG SPACERANGER_SHA256="2566b24f29829b39f3add112a674990b1c54ae2fbe7ccb50a4c7dce9ccf152e6"
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
