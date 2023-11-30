#!/usr/bin/env cached-nix-shell
#! nix-shell -i bash -p fd


# Find modules that have a tests directory
tested=$(fd main.nf.test modules/)

for module in $tested; do
    clean=$(dirname $module | sed 's|tests||')
    rmdir "tests/${clean}"
done
