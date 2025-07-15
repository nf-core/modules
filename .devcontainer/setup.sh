#!/usr/bin/env bash

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Update Nextflow
nextflow self-update

# Install nf-core tools
python -m pip install nf-core

# Install pre-commit hooks
pre-commit install --install-hooks

# Install nf-test
wget -qO- https://get.nf-test.com | bash
sudo mv nf-test /usr/local/bin/nf-test
