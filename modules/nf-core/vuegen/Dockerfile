# Set the base image to the official Python 3.12 slim debian bookworm image
FROM python:3.12-slim-bookworm

# Disable manual feedback for apt-get, update package list, install security updates,
# and install packages (procps, and system dependencies for tinytex and kaleido).
# Then clean up the package list and remove the cache.
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
    procps \
    wget \
    perl \
    libglib2.0-0 \
    libgdk-pixbuf2.0-0 \
    libfontconfig1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install vuegen Python package
RUN pip install --no-cache-dir -U vuegen

# Create and switch to non-root user for security reasons
RUN useradd --create-home appuser
USER appuser

# Set the working directory
WORKDIR /home/appuser

# Install TinyTeX for pdf reports
RUN perl -mFile::Find /dev/null && \
    wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh && \
    ~/.TinyTeX/bin/*/tlmgr install \
    koma-script \
    caption

# Set permissions for the appuser home directory
RUN chmod -R 777 /home/appuser

# Set environment variables for Quarto, dataframe_image, and TinyTeX. Also, get tracebacks from C crashes
ENV XDG_RUNTIME_DIR=/tmp \
    XDG_CACHE_HOME=/tmp/quarto_cache_home \
    XDG_DATA_HOME=/tmp/quarto_data_home \
    PATH="${PATH}:/home/appuser/bin" \
    PYTHONFAULTHANDLER=1
