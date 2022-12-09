# Dockerfile to create container with bcl-convert
# Push to nfcore/bclconvert:<VER>

FROM debian:bullseye-slim
LABEL authors="Matthias De Smet <matthias.desmet@ugent.be>" \
    description="Docker image containing bcl-convert"
# Disclaimer: this container is not provided nor supported by Illumina
# 'ps' command is need by some nextflow executions to collect system stats
# Install procps and clean apt cache
RUN apt-get update \
    && apt-get install -y \
    procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*
# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname
COPY bcl-convert /usr/local/bin/bcl-convert
RUN chmod +x /usr/local/bin/bcl-convert
