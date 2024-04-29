#
# First stage: Quarto installation
#
FROM ubuntu:20.04 as quarto
ARG QUARTO_VERSION=1.3.433
ARG TARGETARCH
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    && apt-get clean

RUN mkdir -p /opt/quarto \
    && curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.tar.gz" \
    && tar -zxvf quarto.tar.gz -C /opt/quarto/ --strip-components=1 \
    && rm quarto.tar.gz

#
# Second stage: Conda environment
#
FROM condaforge/mambaforge:24.3.0-0@sha256:0e77174a59f178c134b2df0982b3fcf5f71d98f31e32b9dd4ff76ec333b186f4
COPY --from=quarto /opt/quarto /opt/quarto
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

LABEL \
    authors = "Erik Fasterius" \
    description = "Dockerfile for the quartonotebook nf-core module"
