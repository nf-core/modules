# Dockerfile to create container with bcl2fastq
# Push to nfcore/bcl2fastq:<VER>

FROM debian:bullseye-slim
LABEL authors="Matthias De Smet <matthias.desmet@ugent.be>" \
    description="Docker image containing bcl2fastq"
# Disclaimer: this container is not provided nor supported by Illumina
# 'ps' command is needed by some nextflow executions to collect system stats
# Install procps and clean apt cache
RUN apt-get update \
    && apt-get install -y \
    procps \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*
# Link hostname cmd to fix hardcoded path
RUN ln -s /bin/hostname /usr/bin/hostname

# Add executable to image
COPY bcl2fastq /usr/local/bin/bcl2fastq

# Add external libs to image
COPY css /usr/local/share/css
COPY xsl /usr/local/share/xsl

# Set permission
RUN chmod +x /usr/local/bin/bcl2fastq && chmod -R +rX /usr/local/share/css && chmod -R +rX /usr/local/share/xsl
