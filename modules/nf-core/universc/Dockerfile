# Dockerfile to create container with UniverSC and Cell Ranger

# DISCLAIMER:
# Cell Ranger is proprietary software owned by 10x Genomics.
# Redistribution is generally restricted by the applicable 10x Genomics
# license and terms. Do not redistribute this image containing Cell Ranger
# unless you are explicitly authorized to do so.
#
# This image is intended for the authorized nf-core use case under the
# applicable agreement with 10x Genomics. This does not grant permission
# to redistribute Cell Ranger independently or outside the scope of that
# authorization.

FROM continuumio/miniconda3:26.7.1-1

ARG CELLRANGER_REDISTRIBUTION_AUTHORIZED=false
RUN if [ "$CELLRANGER_REDISTRIBUTION_AUTHORIZED" != "true" ]; then \
        echo >&2 "ERROR: Building this image bundles Cell Ranger. Confirm redistribution authorization before proceeding."; \
        exit 1; \
    fi

RUN apt-get update \
    && apt-get install -y --no-install-recommends patch procps rename wget \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

RUN conda create -y -p /opt/env -c conda-forge -c bioconda \
        fastq-pair=1.0 \
        multiqc=1.35 \
        samtools=1.16.1 \
    && conda clean -afy

# Pre-downloaded from the 10x Genomics website; ADD extracts it
ENV CELLRANGER_VER=10.1.0
ADD cellranger-${CELLRANGER_VER}.tar.gz /opt/

ENV UNIVERSC_VER=1.2.7
COPY universc.patch /tmp/
RUN wget -qO- https://github.com/minoda-lab/universc/archive/refs/tags/${UNIVERSC_VER}.tar.gz | tar -xz -C /opt \
    && mv /opt/universc-${UNIVERSC_VER} /opt/universc \
    && rm -rf /opt/universc/test \
    && patch -l -d /opt/universc -p1 < /tmp/universc.patch

ENV PATH="/opt/cellranger-${CELLRANGER_VER}:/opt/universc:/opt/env/bin:${PATH}"
