#!/usr/bin/env bash

set -e # fail on error

# include hidden files
# https://askubuntu.com/questions/740805/how-can-i-remove-all-files-from-current-directory-using-terminal
shopt -s dotglob
rm -rf *
