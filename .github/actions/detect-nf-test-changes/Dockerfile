FROM python:3.12-slim-bullseye

COPY requirements.txt requirements.txt

RUN apt-get update \
    && apt-get install -y git \
    && pip install --no-cache-dir --upgrade pip \
    && pip install -r requirements.txt

COPY entrypoint.py /entrypoint.py
COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
