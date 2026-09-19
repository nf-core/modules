# Dockerfile to create container for Vizgen post-processing tool with Cellpose2 plugin
# Push to quay.io/nf-core/vizgen-postprocessing_container:<VER>

FROM python:3.10-slim-bookworm

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLCONFIGDIR=/tmp/mpl_cache \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    CELLPOSE_LOCAL_MODELS_PATH=/tmp/cellpose_models

# Install system libraries
RUN apt-get update && \
    apt-get install -y \
        bzip2 \
        default-jre \
        curl \
        unzip \
        libgl1-mesa-glx \
        libvips42 \
        procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install vpt and plugins
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir vpt vpt-plugin-cellpose2

CMD ["/bin/bash"]
