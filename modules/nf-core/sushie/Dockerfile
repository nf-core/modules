FROM python:3.10

ARG SUSHIE_VERSION=v0.19

RUN apt-get update \
    && apt-get install -y git procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install git+https://github.com/mancusolab/sushie.git@${SUSHIE_VERSION}
