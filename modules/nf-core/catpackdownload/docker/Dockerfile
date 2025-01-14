# Use an official Python runtime as a parent image
FROM python:3.8-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    wget \
    build-essential \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Prodigal
RUN wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux \
    && chmod +x prodigal.linux \
    && mv prodigal.linux /usr/local/bin/prodigal

# Install DIAMOND
RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/diamond \
    && rm diamond-linux64.tar.gz

# Install MMseqs2
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz \
    && tar xvfz mmseqs-linux-sse41.tar.gz \
    && mv mmseqs/bin/mmseqs /usr/local/bin/mmseqs \
    && rm -rf mmseqs*

# Install CAT_pack
WORKDIR /opt
RUN git clone https://github.com/MGXlab/CAT_pack.git

# symlink the tool
RUN ln -sf /opt/CAT_pack/CAT_pack/* /usr/local/bin/
