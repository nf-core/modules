FROM tomkellygenetics/cellranger_clean:3.0.2.9002

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Install system dependencies
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        git-lfs \
        gzip \
        make \
        nano \
        pigz \
        python3-pip \
        rename \
        wget \
        zsh \
    && rm -rf /var/lib/apt/lists/*

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

# Allow Cell Ranger files to be modified
RUN chmod -R 777 /cellranger-3.0.2.9001

# Clone UniverSC, including Git LFS files
RUN git lfs install --system \
    && git clone https://github.com/minoda-lab/universc.git /universc \
    && cd /universc \
    && git lfs pull

# The launcher is now supplied by the cloned repository:
# /universc/launch_universc.sh
# Build the small Cell Ranger reference datasets
RUN cd /universc/test/cellranger_reference/cellranger-tiny-ref \
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
RUN cd /universc \
    && make reference

# Link test FASTQ directories into the Cell Ranger installation
RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-fastq \
    && ln -s \
        /universc/test/shared/cellranger-tiny-fastq/1.2.0 \
        /cellranger-3.0.2.9001/cellranger-tiny-fastq/1.2.0 \
    && ln -s \
        /universc/test/shared/cellranger-tiny-fastq/3.0.0 \
        /cellranger-3.0.2.9001/cellranger-tiny-fastq/3.0.0

# Link test reference directories into the Cell Ranger installation
RUN mkdir -p /cellranger-3.0.2.9001/cellranger-tiny-ref \
    && ln -s \
        /universc/test/cellranger_reference/cellranger-tiny-ref/1.2.0 \
        /cellranger-3.0.2.9001/cellranger-tiny-ref/1.2.0 \
    && ln -s \
        /universc/test/cellranger_reference/cellranger-tiny-ref/3.0.0 \
        /cellranger-3.0.2.9001/cellranger-tiny-ref/3.0.0

# Install the UniverSC launcher in the Cell Ranger executable directory
RUN ln -s \
        /universc/launch_universc.sh \
        /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/bin/universc \
    && ln -s \
        /universc/sub \
        /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/bin/sub \
    && ln -s \
        /universc/whitelists \
        /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/bin/whitelists

# Preserve the original UniverSC compatibility workaround
RUN cp \
        /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/chemistry.py \
        /cellranger-3.0.2.9001/cellranger-cs/3.0.2.9001/lib/python/cellranger/check.py

# Match the permissions used by the original image
RUN chmod -R 777 \
        /cellranger-3.0.2.9001 \
        /universc
