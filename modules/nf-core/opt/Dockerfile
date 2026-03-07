FROM python:3.11-slim

LABEL \
    title="Docker container for package opt - off-target-probe-tracker" \
    maintainer="Sameesh Kher <khersameesh24@gmail.com>"

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/root/.local/bin:$PATH"

# Install build & runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        build-essential \
        pkg-config \
        wget \
        curl \
        procps \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        swig \
        ruby-dev \
        libperl-dev \
        yaggo \
        gettext \
        autoconf \
        automake \
        libtool \
        gffread \
        bowtie2 \
        samtools \
    && rm -rf /var/lib/apt/lists/*

# --- Build and install MUMmer4 v4.0.1 ---
WORKDIR /opt
RUN git clone --branch v4.0.1 https://github.com/mummer4/mummer.git
WORKDIR /opt/mummer
RUN autoreconf -fi && \
    ./configure && \
    make && \
    make install && \
    ldconfig /usr/local/lib
WORKDIR /opt
RUN rm -rf mummer

# --- Install off-target-probe-tracker ---
RUN git clone https://github.com/JEFworks/off-target-probe-tracker.git
WORKDIR /opt/off-target-probe-tracker
RUN pip install --no-cache-dir .
WORKDIR /opt
RUN rm -rf off-target-probe-tracker

# Default work directory
WORKDIR /data
