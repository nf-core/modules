# Dockerfile to create a container with MuSE v

FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git g++ cmake autoconf libtool liblzma-dev zlib1g-dev libbz2-dev libcurl3-dev libssl-dev \
    ca-certificates cpp make libltdl-dev wget unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN git clone --recursive https://github.com/wwylab/MuSE
RUN cd MuSE && ./install_muse.sh

RUN mkdir /MuSE/bin
RUN cp /MuSE/MuSE /MuSE/bin
RUN PATH=$PATH:/MuSE/bin
