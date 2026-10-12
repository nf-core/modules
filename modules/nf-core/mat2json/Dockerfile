FROM ghcr.io/qbic-pipelines/matlabruntime-base:1.0.0@sha256:b1bacc5939d75034a1e9d01f3a294668f5f5b726977d15a7a4beca3e60b7ad0b

ARG MAT2JSON_VERSION="1.0.0"
ARG MAT2JSON_URL="https://github.com/qbic-pipelines/mat2json/releases/download/v1.0.0/mat2json"
ARG MAT2JSON_SHA256="216cac1716c761ebcdde415726762ec462cc79eda7daf3e320549b8765f097e5"

LABEL org.opencontainers.image.authors="Carolin Schwitalla <carolin.schwitalla@qbic.uni-tuebingen.de>"
LABEL org.opencontainers.image.description="Docker image for mat2json, a tool to convert .mat files to .json format, built on MATLAB Runtime R2023a Update 5."
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="$MAT2JSON_VERSION"

RUN mkdir -p /usr/bin/mlrtapp
RUN apt -y update && apt -y install curl
RUN curl -fsSL "$MAT2JSON_URL" -o /usr/bin/mlrtapp/mat2json \
    && echo "$MAT2JSON_SHA256  /usr/bin/mlrtapp/mat2json" | sha256sum -c

RUN chmod -R a+rx /usr/bin/mlrtapp/

# Set the MCR_CACHE_ROOT environment variable
ENV MCR_CACHE_ROOT=/tmp

# Add MATLAB app directory to PATH
ENV PATH="/usr/bin/mlrtapp:${PATH}"
