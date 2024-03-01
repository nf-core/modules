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

# check file integrity using md5sum or vdb-validate as fallback
MD5SUM=$(
	curl -s 'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?filetype=run&acc=!{id}' |
	grep -o '"md5": "[^"]*"' | cut -f4 -d'"'
)
if [ -n "$MD5SUM" ]; then
	[ "$(md5sum !{id}/* | cut -f1 -d' ')" = "$MD5SUM" ]
else
	vdb-validate !{id}
fi || exit 1

cat <<-END_VERSIONS > versions.yml
"!{task.process}":
    sratools: $(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
END_VERSIONS
