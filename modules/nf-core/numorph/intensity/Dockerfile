FROM ghcr.io/qbic-pipelines/matlabruntime-numorphtools:1.0.0@sha256:e05a941552a1f9c741530d68a0d3877363c8e05bb385190f0acade97db73656b

ARG NUMORPH_PREPROCESSING_VERSION="1.0.0"
ARG NUMORPH_PREPROCESSING_URL="https://github.com/qbic-pipelines/Numorph-toolkit/releases/download/1.0.0/numorph_preprocessing"
ARG NUMORPH_PREPROCESSING_SHA256="2adecf1af4a7d362eca7ab2876e017034c4b155df10fa3adcca409a4ac47c976"

LABEL org.opencontainers.image.authors="Carolin Schwitalla <carolin.schwitalla@qbic.uni-tuebingen.de>"
LABEL org.opencontainers.image.description="Docker image for numorph_preprocessing, a collection of tools to preprocess light-sheet microscopy data."
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="$NUMORPH_PREPROCESSING_VERSION"

RUN mkdir -p /usr/bin/mlrtapp
RUN apt -y update && apt -y install curl
RUN curl -fsSL "$NUMORPH_PREPROCESSING_URL" -o /usr/bin/mlrtapp/numorph_preprocessing \
    && echo "$NUMORPH_PREPROCESSING_SHA256  /usr/bin/mlrtapp/numorph_preprocessing" | sha256sum -c

RUN chmod -R a+rx /usr/bin/mlrtapp/

# Set the MCR_CACHE_ROOT environment variable
ENV MCR_CACHE_ROOT=/tmp

# Add MATLAB app directory to PATH
ENV PATH="/usr/bin/mlrtapp:${PATH}"
