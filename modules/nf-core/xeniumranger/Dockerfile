#  Use a base image with the required operating system
FROM ubuntu:20.04 AS builder

LABEL authors="Sameesh Kher <sameesh.kher@dkfz-heidelberg.de>" \
    description="Docker image containing Xenium Ranger"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

ARG VERSION="2.0.0"

# copy over xeniumranger
# xenium ranger has been downloaded from https://www.10xgenomics.com/support/software/xenium-ranger/downloads
COPY xeniumranger-${VERSION}.tar /opt/xeniumranger-${VERSION}.tar

# install xenium ranger
RUN tar -xzvf /opt/xeniumranger-${VERSION}.tar -C /opt/ && \
    rm /opt/xeniumranger-${VERSION}.tar

# Set environment variables
ENV PATH="/opt/xeniumranger-xenium2.0/bin:$PATH"

# multistage to reduce image size
FROM ubuntu:20.04

# Set environment variables
ENV PATH="/opt/xeniumranger-xenium2.0/bin:$PATH"

# copy over xenium from builder
COPY --from=builder /opt/xeniumranger-xenium2.0/ /opt/xeniumranger-xenium2.0/