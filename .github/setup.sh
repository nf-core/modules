#!/bin/bash

## Setup script from https://github.com/bioconda/bioconda-recipes/blob/e5e9477c8e24e6edaa915e0042bc3531d04b2ec7/.circleci/setup.sh
set -e
[[ -z $WORKSPACE ]] && WORKSPACE=`pwd`
[[ -z $BOOTSTRAP ]] && BOOTSTRAP=false
[[ -z $BASH_ENV ]] && BASH_ENV=`mktemp`
[[ -z $USE_DOCKER ]] && USE_DOCKER=true

set -u

source .github/common.sh

cat >> $BASH_ENV <<EOF

# Set path
export PATH="${WORKSPACE}/miniconda/bin:${PATH}"

if [ -f "${WORKSPACE}/miniconda/etc/profile.d/conda.sh" ] ; then
    . "${WORKSPACE}/miniconda/etc/profile.d/conda.sh"
fi
EOF

. $BASH_ENV

if ! type bioconda-utils 2> /dev/null || [[ $BOOTSTRAP == "true" ]]; then
    echo "Setting up bioconda-utils..."

    # setup conda and bioconda-utils if not loaded from cache
    mkdir -p $WORKSPACE

    # step 1: download and install miniconda
    if [[ $OSTYPE == darwin* ]]; then
        tag="MacOSX"
    elif [[ $OSTYPE == linux* ]]; then
        tag="Linux"
    else
        echo "Unsupported OS: $OSTYPE"
        exit 1
    fi
    curl -L -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-$MINICONDA_VER-$tag-x86_64.sh
    bash miniconda.sh -b -p $WORKSPACE/miniconda
    . $BASH_ENV
    conda activate base

    # step 2: setup channels
    conda config --system --add channels defaults
    conda config --system --add channels bioconda
    conda config --system --add channels conda-forge

    # step 3: install bioconda-utils
    additional_packages='git pip'
    if [[ $OSTYPE == darwin* ]]; then
        # Pinned to 2.5.3 to make sure we don't get unexpected changes.
        additional_packages="${additional_packages} conda-forge-ci-setup=2.5.3"
    fi
    conda install -y $additional_packages --file https://raw.githubusercontent.com/bioconda/bioconda-utils/$BIOCONDA_UTILS_TAG/bioconda_utils/bioconda_utils-requirements.txt
    pip install git+https://github.com/bioconda/bioconda-utils.git@$BIOCONDA_UTILS_TAG

    # step 4: configure local channel
    mkdir -p $WORKSPACE/miniconda/conda-bld/{noarch,linux-64,osx-64}
    conda index $WORKSPACE/miniconda/conda-bld
    conda config --system --add channels file://$WORKSPACE/miniconda/conda-bld

    # step 5: cleanup
    conda clean -y --all
    rm miniconda.sh
fi

# # Fetch the master branch for comparison (this can fail locally, if git remote 
# # is configured via ssh and this is executed in a container).
# if [[ $BOOTSTRAP != "true" ]]; then
#     git fetch origin +master:master || true
# fi
