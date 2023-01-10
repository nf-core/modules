FROM ubuntu:22.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        ffmpeg \
        git \
        openjdk-11-jdk-headless \
        python3-dev \
        python3-pip \
        python3-tk \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -q -U \
    numpy \
    pip

COPY / /app/ashlar/
RUN pip3 install /app/ashlar

ENV OMP_NUM_THREADS 1

VOLUME /data
VOLUME /data2
VOLUME /data3

WORKDIR /data
