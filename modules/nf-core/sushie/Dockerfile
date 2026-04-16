FROM python:3.10-slim

ARG SUSHIE_VERSION=v0.19

RUN apt-get update \
    && apt-get install -y gcc git procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir git+https://github.com/mancusolab/sushie.git@${SUSHIE_VERSION}
