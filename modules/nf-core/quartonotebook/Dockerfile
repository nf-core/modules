FROM condaforge/mambaforge:24.9.2-0
#@sha256:e4edd160d0a10c6f577b127021d853b460e7e8c1bf1850c9aa4a52c4b8d9c8ee
ARG QUARTO_VERSION=1.3.433
ARG TARGETARCH

LABEL org.opencontainers.image.authors="Erik Fasterius"
LABEL org.opencontainers.image.description="Dockerfile for the quartonotebook nf-core module"
LABEL org.opencontainers.image.version="$QUARTO_VERSION"
LABEL org.opencontainers.image.documentation="https://github.com/nf-core/modules/blob/master/modules/nf-core/quartonotebook/README.md"
LABEL org.opencontainers.image.source="https://github.com/nf-core/modules"
LABEL org.opencontainers.image.vendor="nf-core"
LABEL org.opencontainers.image.license="https://github.com/nf-core/modules/blob/master/LICENSE"

ADD https://github.com/quarto-dev/quarto-cli#v${QUARTO_VERSION} /opt/quarto
ENV PATH="${PATH}:/opt/quarto/bin"

# Install packages using Mamba; also remove static libraries, python bytecode
# files and javascript source maps that are not required for execution
COPY environment.yml ./
RUN mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && find /opt/conda -follow -type f -name '*.js.map' -delete

CMD /bin/bash
