#!/usr/bin/env bash

set -u


# Get the expected NCBI settings path and define the environment variable
# `NCBI_SETTINGS`.
eval "\$(vdb-config -o n NCBI_SETTINGS | sed 's/[" ]//g')"

# If the user settings do not exist yet, create a file suitable for `prefetch`
# and `fasterq-dump`. Otherwise, copy the existing settings file.
if [[ ! -f "\${NCBI_SETTINGS}" ]]; then
    printf "$config" > 'user-settings.mkfg'
else
    if [[ "\${NCBI_SETTINGS}" != *.mkfg ]]; then
        echo "The detected settings '\${NCBI_SETTINGS}' do not have the required"\
        "file extension '.mkfg'." >&2
        exit 1
    fi
    cp "\${NCBI_SETTINGS}" ./
fi
