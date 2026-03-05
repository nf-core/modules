FROM matlabruntime/r2023a/release/update5/300000000000000000

ARG MAT2JSON_VERSION="1.0.0"
ARG MAT2JSON_URL="https://github.com/qbic-pipelines/mat2json/releases/download/v1.0.0/mat2json"
ARG MAT2JSON_SHA256="216cac1716c761ebcdde415726762ec462cc79eda7daf3e320549b8765f097e5"

LABEL org.opencontainers.image.authors="Carolin Schwitalla <carolin.schwitalla@qbic.uni-tuebingen.de>"
LABEL org.opencontainers.image.description="Docker image for mat2json, a tool to convert .mat files to .json format, built on MATLAB Runtime R2023a Update 5."
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="$MAT2JSON_VERSION"

RUN mkdir -p /usr/bin/mlrtapp

ADD --checksum=sha256:$MAT2JSON_SHA256 \
    $MAT2JSON_URL \
    /usr/bin/mlrtapp

RUN chmod -R a+rx /usr/bin/mlrtapp/

# Set the MCR_CACHE_ROOT environment variable
ENV MCR_CACHE_ROOT=/tmp

# Add MATLAB app directory to PATH
ENV PATH="/usr/bin/mlrtapp:${PATH}"
