#!/usr/bin/env bash

# Customise the terminal command prompt
echo "export PROMPT_DIRTRIM=2" >> $HOME/.bashrc
echo "export PS1='\[\e[3;36m\]\w ->\[\e[0m\\] '" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]\w ->\[\e[0m\\] '

# Update Nextflow
nextflow self-update

# Install prek hooks
pip install prek
prek install --install-hooks

# Update welcome message
echo "Welcome to nf-core/modules devcontainer!" > /usr/local/etc/vscode-dev-containers/first-run-notice.txt
