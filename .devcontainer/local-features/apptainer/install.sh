#!/usr/bin/env bash

# Install Apptainer (Singularity)
add-apt-repository -y ppa:apptainer/ppa
apt-get update --quiet

apt install -y apptainer

apt-get clean
rm -rf /var/lib/apt/lists/*
