FROM nfcore/base:2.1
# Install the conda environment


COPY ./somatosim/environment.yml /
ENV CONDA_VERBOSITY=2

RUN conda env create -f /environment.yml -v && conda clean -a
ENV PATH /opt/conda/envs/nf-core-somatosim/bin:$PATH

# Add conda installation dir to PATH (instead of doing 'conda activate')


RUN conda install -c bioconda "samtools>=1.10"
RUN conda install pip

RUN git clone https://github.com/BieseckerLab/SomatoSim.git

