FROM ghcr.io/qbic-pipelines/matlabruntime-numorphtools:1.0.0@sha256:e05a941552a1f9c741530d68a0d3877363c8e05bb385190f0acade97db73656b

ARG NUMORPH_ANALYZE_VERSION="1.0.0"
ARG NUMORPH_ANALYZE_URL="https://github.com/qbic-pipelines/Numorph-toolkit/releases/download/1.0.0/numorph_analyze"
ARG NUMORPH_ANALYZE_SHA256="59e88b81ff0369aa00d982699641b86a9874b2b4f59d403b64bbe8127b1d25d1"

LABEL org.opencontainers.image.authors="Carolin Schwitalla <carolin.schwitalla@qbic.uni-tuebingen.de>"
LABEL org.opencontainers.image.description="Docker image for numorph_analyze, a collection of tools to register light-sheet microscopy images to the ara reference atlas."
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="$NUMORPH_ANALYZE_VERSION"

RUN mkdir -p /usr/bin/mlrtapp
RUN apt -y update && apt -y install curl
RUN curl -fsSL "$NUMORPH_ANALYZE_URL" -o /usr/bin/mlrtapp/numorph_analyze \
    && echo "$NUMORPH_ANALYZE_SHA256  /usr/bin/mlrtapp/numorph_analyze" | sha256sum -c

RUN chmod -R a+rx /usr/bin/mlrtapp/

# Set the MCR_CACHE_ROOT environment variable
ENV MCR_CACHE_ROOT=/tmp

# Add MATLAB app directory to PATH
ENV PATH="/usr/bin/mlrtapp:${PATH}"
