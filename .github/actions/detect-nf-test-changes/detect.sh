#!/usr/bin/env bash
# Usage: detect.sh <base-ref> [include-tags] [exclude-tags]
# Prints a JSON array of nf-core component directories whose nf-tests are related
# to the changes between merge-base(<base-ref>, HEAD) and HEAD.
set -euo pipefail

base_ref="$1"
include_tags="${2:-}"
exclude_tags="${3:-}"

base=$(git merge-base "$base_ref" HEAD)

tests=$(nf-test test --dry-run --changed-since "$base" --related-tests --filter process,workflow 2>&1 \
    | awk '/related test\(s\)/ {f=1; next} /Dry run mode/ {f=0} f && /\.nf\.test$/ {sub(/^[[:space:]]*•[[:space:]]*/, ""); print}' \
    | sed "s#^$PWD/##")

has_tag() { grep -qE "^[[:space:]]*tag[[:space:]]+[\"']$2[\"']" "$1"; }

components=()
for t in $tests; do
    keep=1
    if [ -n "$include_tags" ]; then
        keep=0
        for tag in ${include_tags//,/ }; do has_tag "$t" "$tag" && keep=1; done
    fi
    for tag in ${exclude_tags//,/ }; do has_tag "$t" "$tag" && keep=0; done
    [ "$keep" = 1 ] || continue
    comp=$(dirname "$(dirname "$t")")
    [[ "$comp" =~ ^(modules|subworkflows)/nf-core/ ]] && components+=("$comp")
done

printf '%s\n' "${components[@]+"${components[@]}"}" | { grep . || true; } | sort -u | jq -R . | jq -cs .
