#!/usr/bin/env bash
# Bump tool versions and move nf-core modules from quay.io to Seqera Containers.
# Usage: convert.sh <module> [<module> ...]   (run from the nf-core/modules root)
# Runs sequentially: parallel nf-core runs race on ~/.config/nfcore's git cache.
set -uo pipefail

log_dir="${TMPDIR:-/tmp}/seqera-containers"
mkdir -p "$log_dir"
failed=()

for m in "$@"; do
    dir="modules/nf-core/$m"
    log="$log_dir/$(echo "$m" | tr / _).log"
    [[ -f "$dir/main.nf" ]] || { echo "SKIP $m: no $dir/main.nf"; failed+=("$m"); continue; }

    nf-core modules bump-versions "$m" > "$log" 2>&1
    # bump-versions builds containers itself when it bumps; only create when still on quay.
    if grep -qE "quay.io|depot.galaxyproject.org" "$dir/main.nf"; then
        nf-core modules containers create "$m" >> "$log" 2>&1
    fi

    if grep -qE "quay.io|depot.galaxyproject.org" "$dir/main.nf"; then
        echo "FAIL $m (see $log)"
        failed+=("$m")
    else
        echo "OK   $m $(grep -o "community.wave.seqera.io/library/[^'\"]*" "$dir/main.nf")"
    fi
done

if (( ${#failed[@]} )); then
    echo "Failed: ${failed[*]}"
    exit 1
fi
