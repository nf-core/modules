#  Use a base image with the required operating system
FROM ubuntu:20.04 AS builder

LABEL authors="Sameesh Kher <sameesh.kher@dkfz-heidelberg.de>" \
    description="Docker image containing Xenium Ranger"

RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# copy over xeniumranger
# the latest version of the xeniumranger tool has been downloaded from https://www.10xgenomics.com/support/software/xenium-ranger/downloads
COPY xeniumranger-4.0.0.tar.gz /xeniumranger-4.0.0.tar.gz

# install xenium ranger
RUN tar -xzvf /xeniumranger-4.0.0.tar.gz && \
    rm /xeniumranger-4.0.0.tar.gz

# Set environment variables
# ENV PATH="xeniumranger-xenium4.0/bin:$PATH"

# multistage to reduce image size
FROM ubuntu:20.04

# Set environment variables
ENV PATH="/xeniumranger-xenium4.0/bin:$PATH"

# copy over xenium from builder
COPY --from=builder /xeniumranger-xenium4.0 /xeniumranger-xenium4.0/
