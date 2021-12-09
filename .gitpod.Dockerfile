FROM gitpod/workspace-full as builder

USER root

# Install custom tools, runtime, etc.
RUN apt-get update && apt-get install -y \
    build-essential \
    uuid-dev \
    libgpgme-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup-bin && rm -rf /var/lib/apt/lists/*

# Install Singularity (Go is already installed)
RUN wget https://github.com/sylabs/singularity/releases/download/v3.9.1/singularity-ce-3.9.1.tar.gz && \
    tar -xzf singularity-ce-3.9.1.tar.gz && \
    cd singularity-ce-3.9.1 && \
    ./mconfig --without-suid -p /usr/local/singularity && \
    make -C ./builddir && \
    make -C ./builddir install


FROM gitpod/workspace-full

# Copy over singularity and setup dependencies
COPY --from=builder /usr/local/singularity /usr/local/singularity
ENV PATH="/usr/local/singularity/bin:$PATH" \
    SINGULARITY_TMPDIR="/tmp-singularity"
RUN sudo apt-get update && sudo apt-get install -y \
    squashfs-tools && \
    sudo rm -rf /var/lib/apt/lists/* && \
    sudo mkdir -p $SINGULARITY_TMPDIR

# Install Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sudo bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN sudo chown -R gitpod:gitpod /opt/conda $SINGULARITY_TMPDIR

# Install Mamba
RUN conda install mamba -n base -c bioconda -c conda-forge && \
    conda clean --all -f -y

# Install Nextflow nf-core pytest-workflow
RUN mamba install nextflow=21.10.0 nf-core=2.1 pytest-workflow=1.6.0 -n base -c bioconda -c conda-forge && \
    mamba clean --all -f -y

