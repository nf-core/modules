# ---- Base image ----
FROM python:3.14

LABEL authors="Clara Bierwirth" \
      description="Stitchr container with IMGT human and mouse data"

# ---- Environment settings ----
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/usr/local/bin:${PATH}"

# ---- System dependencies ----
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      build-essential \
      wget \
    && rm -rf /var/lib/apt/lists/*

# ---- Install Stitchr and IMGTgeneDL ----
RUN pip install --no-cache-dir --upgrade --force-reinstall \
    stitchr==1.3.1 \
    IMGTgeneDL==0.6.1

# ---- Download human and mouse IMGT data ----
RUN stitchrdl -s human
RUN stitchrdl -s mouse
