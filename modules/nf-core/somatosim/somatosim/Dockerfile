FROM ubuntu:20:04

RUN apt -y update
RUN apt -y install git

RUN git clone https://github.com/BieseckerLab/SomatoSim.git
RUN cd SomatoSim && python3 -m pip install . 

