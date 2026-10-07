# Dockerfile to create container with Cell Ranger
# Push to quay.io/nf-core/universc

FROM continuumio/miniconda3:latest

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install system dependencies
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        git \
        git-lfs \
        gzip \
        pigz \
        rename \
        cpio \
        procps \
        rpm2cpio \
        unzip \
        wget \
        curl \
        zsh \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Set different version
ENV CELLRANGER_VER=10.1.0
ENV CELLRANGER_HOME=/opt/cellranger-${CELLRANGER_VER}
ENV UNIVERSC_HOME=/opt/universc

# TODO Update CELLRANGE_VER and URL cellranger file
WORKDIR /opt
RUN wget -O cellranger-$CELLRANGER_VER.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VER}.tar.gz?Expires=1791420039&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=iKFpqhsq5MZZ6S~dJX0M8xi1LjO7zfcV-7tg1Hm2PL4nXE~KMf-oT15B6gkLIagps8DTditel79Gz~P09PBU~e3LStXs3bBKgloZGKwtsQJ2oSgvMaCz7kz0r1qycW2H2sgLIEC7dFW6eA3WOvaN86-llRE9A4HHFfLGlJJB7niIsKPo34VeQ8KaMmO3q3AUqyyF4QurDK3Z6zHhxSyvFhGDCDaZ~gnft-~UKk~1rS2g2P0TCc6m-CTthBlom6mgk5PBzV~lOmEaTpIjEDSGbuSyaLs4Yi1a5jT6StC2RyKqDAM0mZcViS3Iu2FIq3D0iFMlq6PxbJ61UvMA2kUw8Q__" \
    && tar -xzf "cellranger-$CELLRANGER_VER.tar.gz" \
    && test -x "$CELLRANGER_HOME/cellranger" \
    && rm -f "cellranger-$CELLRANGER_VER.tar.gz"

ENV PATH="${CELLRANGER_HOME}:${CELLRANGER_HOME}/bin:${UNIVERSC_HOME}:${PATH}"

# Install Oh My Zsh non-interactively
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" \
        "" --unattended

# Install MultiQC
RUN python3 -m pip install --no-cache-dir \
        --upgrade pip setuptools wheel \
    && python3 -m pip install --no-cache-dir multiqc

# Install fastq-pair
RUN git clone https://github.com/linsalrob/fastq-pair.git /tmp/fastq-pair \
    && cd /tmp/fastq-pair \
    && mkdir build \
    && cd build \
    && gcc -std=gnu99 \
        ../main.c \
        ../robstr.c \
        ../fastq_pair.c \
        ../is_gzipped.c \
        -o fastq_pair \
    && install -m 0755 fastq_pair /usr/local/bin/fastq_pair \
    && rm -rf /tmp/fastq-pair

# Install STAR aligner
RUN wget -q \
        https://github.com/alexdobin/STAR/archive/refs/tags/2.5.1b.tar.gz \
        -O /tmp/STAR-2.5.1b.tar.gz \
    && tar -xzf /tmp/STAR-2.5.1b.tar.gz -C /tmp \
    && cd /tmp/STAR-2.5.1b \
    && make \
    && install -m 0755 bin/Linux_x86_64/STAR* /usr/local/bin/ \
    && cd source \
    && make \
    && install -m 0755 STAR /usr/local/bin/STAR \
    && rm -rf /tmp/STAR-2.5.1b /tmp/STAR-2.5.1b.tar.gz

# Clone UniverSC, including Git LFS files
RUN git lfs install --system \
    && git clone https://github.com/minoda-lab/universc.git "${UNIVERSC_HOME}" \
    && cd "${UNIVERSC_HOME}" \
    && git lfs pull

# The launcher is now supplied by the cloned repository:
# /universc/launch_universc.sh
# Build the small Cell Ranger reference datasets
RUN cd "${UNIVERSC_HOME}/test/cellranger_reference/cellranger-tiny-ref" \
    && rm -rf 3.0.0 1.2.0 \
    && cellranger mkref \
        --genome=3.0.0 \
        --fasta=genome-3.0.0.fa \
        --genes=genes-3.0.0.gtf \
    && cellranger mkref \
        --genome=1.2.0 \
        --fasta=genome-1.2.0.fa \
        --genes=genes-1.2.0.gtf

# Build UniverSC references
RUN cd ${UNIVERSC_HOME} \
    && make reference

# Link test FASTQ directories into the Cell Ranger installation
RUN ln -s \
        "${UNIVERSC_HOME}/test/shared/cellranger-tiny-fastq" \
        "${CELLRANGER_HOME}/cellranger-tiny-fastq" \
    && ln -s \
        "${UNIVERSC_HOME}/test/cellranger_reference/cellranger-tiny-ref" \
        "${CELLRANGER_HOME}/cellranger-tiny-ref"

# Install the UniverSC launcher in the Cell Ranger executable directory
RUN ln -s \
        "${UNIVERSC_HOME}/launch_universc.sh" \
        "${CELLRANGER_HOME}/bin/universc" \
    && ln -s \
        "${UNIVERSC_HOME}/sub" \
        "${CELLRANGER_HOME}/bin/sub" \
    && ln -s \
        "${UNIVERSC_HOME}/whitelists" \
        "${CELLRANGER_HOME}/bin/whitelists"

# Preserve the original UniverSC compatibility workaround
RUN cp \
        "${CELLRANGER_HOME}/lib/python/cellranger/chemistry.py" \
        "${CELLRANGER_HOME}/lib/python/cellranger/check.py"

# FIX sed -i 's/"\$technology" == "vasa-drop"/ "\$technology" == "vasa-drop"/' "\$local_universc/universc"

# Match the permissions used by the original image
RUN chmod -R a+rwX \
    "${CELLRANGER_HOME}" \
    "${UNIVERSC_HOME}"
