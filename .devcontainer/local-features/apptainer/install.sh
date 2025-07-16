#!/usr/bin/env bash

# Install Apptainer (Singularity)

apt-get update --quiet

# installs add-apt-repository
apt install --reinstall -y software-properties-common

add-apt-repository -y ppa:apptainer/ppa

apt install -y apptainer

apt-get clean
rm -rf /var/lib/apt/lists/*
