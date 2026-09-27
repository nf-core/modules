#!/usr/bin/env bash
# Usage: detect.sh <base-ref> [include-tags] [exclude-tags]
# Prints a JSON array of nf-test files related to the changes between merge-base(<base-ref>, HEAD) and HEAD.
set -euo pipefail

base_ref="$1"
include_tags="${2:-}"
exclude_tags="${3:-}"

base=$(git merge-base "$base_ref" HEAD)
csv=$(mktemp)
rm -f "$csv"

log=$(nf-test test --dry-run --changed-since "$base" --related-tests --filter process,workflow --csv "$csv" 2>&1) || {
    echo "$log" >&2
    echo "nf-test dry run failed" >&2
    exit 1
}

if [ -s "$csv" ]; then
    tests=$(tail -n +2 "$csv" | cut -d'"' -f2 | sed "s#^$PWD/##" | sort -u)
elif grep -q "Nothing to do" <<<"$log"; then
    tests=""
else
    echo "$log" >&2
    echo "unrecognised nf-test output: neither a CSV report nor 'Nothing to do'" >&2
    exit 1
fi
rm -f "$csv"

has_tag() { grep -qE "^[[:space:]]*tag[[:space:]]+[\"']$2[\"']" "$1"; }

selected=()
for t in $tests; do
    [[ "$t" =~ ^(modules|subworkflows)/nf-core/ ]] || continue
    keep=1
    if [ -n "$include_tags" ]; then
        keep=0
        for tag in ${include_tags//,/ }; do has_tag "$t" "$tag" && keep=1; done
    fi
    for tag in ${exclude_tags//,/ }; do has_tag "$t" "$tag" && keep=0; done
    [ "$keep" = 1 ] && selected+=("$t")
done

printf '%s\n' "${selected[@]+"${selected[@]}"}" | { grep . || true; } | jq -R . | jq -cs .
