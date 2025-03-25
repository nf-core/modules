# Use an official Python image as the base
FROM python:3.11-slim

# Install system dependencies (including `procps` for `ps` command)
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Set a working directory inside the container
WORKDIR /app

# Copy the requirements file into the container
COPY requirements.txt .

# Install Python packages
RUN pip install --no-cache-dir -r requirements.txt

# Ensure the container runs as a non-root user for security
RUN groupadd -r appuser && useradd --no-log-init -r -g appuser appuser
USER appuser
