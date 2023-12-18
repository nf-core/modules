#!/usr/bin/env cached-nix-shell
#! nix-shell -i bash -p fd yq-go


# Find modules that have a tests directory
tested=$(fd main.nf.test modules/)

for module in $tested; do
    clean=$(dirname $module | sed 's|/tests||' | sed 's|modules/nf-core/||')
    yq -i "del(.${clean})" tests/config/pytest_modules.yml
    # rm -rf "tests/${clean}"
done
