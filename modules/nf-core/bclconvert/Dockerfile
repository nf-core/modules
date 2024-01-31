# Dockerfile to create container with bcl-convert
# Push to nfcore/bclconvert:<VER>

FROM debian:bullseye-slim
LABEL authors="Matthias De Smet <matthias.desmet@ugent.be>" \
    description="Docker image containing bcl-convert"
ARG BCL_VERSION="4.2.7"
# Disclaimer: this container is not provided nor supported by Illumina
# 'ps' command is need by some nextflow executions to collect system stats
# Install procps and clean apt cache
RUN apt-get update \
    && apt-get install -y \
    procps \
    rpm2cpio \
    cpio \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*
# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname
# Install bcl-convert
ADD https://s3.amazonaws.com/webdata.illumina.com/downloads/software/bcl-convert/bcl-convert-${BCL_VERSION}-2.el7.x86_64.rpm bcl-convert.rpm
RUN rpm2cpio bcl-convert.rpm | cpio -i --make-directories
