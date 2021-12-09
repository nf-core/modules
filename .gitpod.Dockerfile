# Singularity installation adopted from https://github.com/tyson-swetnam/singularity-gitpod/blob/main/.gitpod.Dockerfile
FROM gitpod/workspace-full

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

# create /workspace in container (this is also done by k8s at initiation), change owner to chown
# RUN mkdir /workspace && chown -R 33333:33333 /workspace

# Install Singularity (Go is already installed)
RUN wget https://github.com/sylabs/singularity/releases/download/v3.9.1/singularity-ce-3.9.1.tar.gz && \
    tar -xzf singularity-ce-3.9.1.tar.gz && \
    cd singularity-ce-3.9.1 && \
    ./mconfig && \
    make -C ./builddir && \
    make -C ./builddir install && \
    cd - && rm -rf singularity-ce-3.9.* && \   
    echo ". /usr/local/etc/bash_completion.d/singularity" >> ${HOME}/.bashrc

# Install Conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN chown -R gitpod:gitpod /opt/conda 

USER gitpod

# Install nf-core, Mamba, and pytest-workflow
RUN conda install nextflow nf-core pytest-workflow mamba -n base -c conda-forge -c bioconda && \
    conda clean --all -f -y


