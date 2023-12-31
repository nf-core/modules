FROM nfcore/base:2.1

# Install the conda environment

# COPY environment.yml /
# RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
RUN conda install -c bioconda "samtools>=1.18"


RUN git clone https://github.com/BieseckerLab/SomatoSim.git

RUN cd SomatoSim && python3 -m pip install . 

