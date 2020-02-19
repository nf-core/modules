#!/bin/bash
WORKSPACE=$(pwd)

## Common definitions for the bioconda build system. 
MINICONDA_VER=4.6.14
BIOCONDA_UTILS_TAG=master

if ! type bioconda-utils 2> /dev/null || [[ $BOOTSTRAP == "true" ]]; then
    echo "Setting up bioconda-utils..."

    # setup channels
    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge

    # install bioconda-utils
    # fix temporary problem with bioconda-utlis not being compatible with the 
    # latest version of gitdb
    additional_packages="gitdb2<3"
    conda install -y $additional_packages --file https://raw.githubusercontent.com/grst/bioconda-utils/$BIOCONDA_UTILS_TAG/bioconda_utils/bioconda_utils-requirements.txt
    pip install git+https://github.com/grst/bioconda-utils.git@$BIOCONDA_UTILS_TAG

    # step 4: configure local channel
    mkdir -p $WORKSPACE/miniconda/conda-bld/{noarch,linux-64,osx-64}
    conda index $WORKSPACE/miniconda/conda-bld
    conda config --system --add channels file://$WORKSPACE/miniconda/conda-bld

    # step 5: cleanup
    conda clean -y --all
fi

# # Fetch the master branch for comparison (this can fail locally, if git remote 
# # is configured via ssh and this is executed in a container).
# if [[ $BOOTSTRAP != "true" ]]; then
#     git fetch origin +master:master || true
# fi
