#!/usr/bin/env bash

set -u

retry_with_backoff() {
    local max_attempts=${1}
    local delay=${2}
    local max_time=${3}
    local attempt=1
    local output=
    local status=

    # Remove the first three arguments to this function in order to access
    # the 'real' command with `${@}`.
    shift 3

    while [ ${attempt} -le ${max_attempts} ]; do
        output=$("${@}")
        status=${?}

        if [ ${status} -eq 0 ]; then
            break
        fi

        if [ ${attempt} -lt ${max_attempts} ]; then
            echo "Failed attempt ${attempt} of ${max_attempts}. Retrying in ${delay} s." >&2
            sleep ${delay}
        elif [ ${attempt} -eq ${max_attempts} ]; then
            echo "Failed after ${attempt} attempts." >&2
            return ${status}
        fi

        attempt=$(( ${attempt} + 1 ))
        delay=$(( ${delay} * 2 ))
        if [ ${delay} -ge ${max_time} ]; then
            delay=${max_time}
        fi
    done

    echo "${output}"
}

export NCBI_SETTINGS="$PWD/!{ncbi_settings}"

retry_with_backoff !{args2} \
    prefetch \
    !{args} \
    !{id}

# check file integrity using vdb-validate or (when archive contains no checksums) md5sum
vdb-validate !{id} > vdb-validate_result.txt 2>&1 || exit 1
if grep -q "checksums missing" vdb-validate_result.txt; then
    VALID_MD5SUMS=$(curl --silent --fail --location --retry 3 --retry-delay 60 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&acc=!{id}')
    LOCAL_MD5SUMS=$(md5sum !{id}/* | cut -f1 -d' ')
    if ! grep -q -F -f <(echo "$LOCAL_MD5SUMS") <(echo "$VALID_MD5SUMS"); then
        echo "MD5 sum check failed" 1>&2
        exit 1
    fi
fi

cat <<-END_VERSIONS > versions.yml
"!{task.process}":
    sratools: $(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    curl: $(curl --version | head -n 1 | sed 's/^curl //; s/ .*$//')
END_VERSIONS
