FROM nfcore/base:2.1
LABEL \
    author="Maxime Garcia" \
    description="snpEff image for nf-core pipelines" \
    maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment

# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
RUN conda install -c bioconda "samtools>=1.10"


RUN git clone https://github.com/BieseckerLab/SomatoSim.git

RUN cd SomatoSim && python3 -m pip install . 

