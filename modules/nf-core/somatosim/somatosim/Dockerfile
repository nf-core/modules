FROM quay.io/biocontainers/samtools:0.1.19--2

RUN apt-get update
RUN apt-get -y update && apt-get -y install git

RUN git clone https://github.com/BieseckerLab/SomatoSim.git
WORKDIR SomatoSim
RUN cd SomatoSim && python3 -m pip install . 

