# Using the Ubuntu image (our OS)
FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
# Update package manager (apt-get) 
# and install (with the yes flag `-y`)
# Python and Pip
RUN apt-get update && apt-get install -y \
    python3.11  \
    python3-pip \
    wget git
RUN ln -sf python3.11 /usr/bin/python3

# Add conda installation dir to PATH (instead of doing 'conda activate')
 


RUN wget "http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda install -c conda-forge libgcc-ng pip
RUN conda install -c bioconda "samtools=1.18"




RUN git clone https://github.com/BieseckerLab/SomatoSim.git

RUN cd SomatoSim && pip3 install .
