#  Use a base image with the required operating system
FROM ubuntu:20.04 AS builder

LABEL authors="Sameesh Kher <sameesh.kher@dkfz-heidelberg.de>" \
    description="Docker image containing Xenium Ranger"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

ARG MAJOR="2"
ARG MINOR="0"
ARG PATCH="0"

# copy over xeniumranger
# the latest version of the xeniumranger tool has been downloaded from https://www.10xgenomics.com/support/software/xenium-ranger/downloads
COPY xeniumranger-${MAJOR}.${MINOR}.${PATCH}.tar /opt/xeniumranger-${MAJOR}.${MINOR}.${PATCH}.tar

# install xenium ranger
RUN tar -xzvf /opt/xeniumranger-${MAJOR}.${MINOR}.${PATCH}.tar -C /opt/ && \
    rm /opt/xeniumranger-${MAJOR}.${MINOR}.${PATCH}.tar

# Set environment variables
ENV PATH="/opt/xeniumranger-xenium${MAJOR}.${MINOR}/bin:$PATH"

# multistage to reduce image size
FROM ubuntu:20.04

# Set environment variables
ENV PATH="/opt/xeniumranger-xenium${MAJOR}.${MINOR}/bin:$PATH"

# copy over xenium from builder
COPY --from=builder /opt/xeniumranger-xenium${MAJOR}.${MINOR}/ /opt/xeniumranger-xenium${MAJOR}.${MINOR}/