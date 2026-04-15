# =============================================
# FlashPCA 2.0 - Dockerfile per nf-core
# Basato sul Makefile ufficiale di Gad Abraham (GPL-3.0)
# Creato per nf-core da dbaku42
# =============================================

FROM ubuntu:22.04

LABEL maintainer="dbaku42"
LABEL description="FlashPCA - Fast Principal Component Analysis of large-scale genomic data (GPL-3.0)"
LABEL org.opencontainers.image.source="https://github.com/gabraham/flashpca"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    ca-certificates \
    libeigen3-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Download Spectra (required dependency)
RUN wget -q https://github.com/yixuan/spectra/archive/refs/tags/v0.8.1.tar.gz -O /tmp/spectra-0.8.1.tar.gz && \
    tar -xzf /tmp/spectra-0.8.1.tar.gz -C /tmp && \
    rm /tmp/spectra-0.8.1.tar.gz

# Clone and build FlashPCA
RUN git clone --depth 1 https://github.com/gabraham/flashpca.git /opt/flashpca && \
    cd /opt/flashpca && \
    make -j$(nproc) all \
        EIGEN_INC=/usr/include/eigen3 \
        BOOST_INC=/usr/include/boost \
        SPECTRA_INC=/tmp/spectra-0.8.1/include && \
    cp flashpca /usr/local/bin/flashpca && \
    chmod +x /usr/local/bin/flashpca && \
    ln -s /usr/local/bin/flashpca /usr/local/bin/flashpca_x86-64

# Cleanup
RUN rm -rf /opt/flashpca /tmp/spectra-0.8.1 /var/cache/* /root/.cache

WORKDIR /data

ENTRYPOINT ["flashpca"]
CMD ["--help"]

# Version info for nf-core
RUN echo "FlashPCA version information:" && \
    /usr/local/bin/flashpca --version 2>&1 || echo "FlashPCA 2.0 (built from source)"
