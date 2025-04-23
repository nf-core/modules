#!/usr/bin/env bash

# Find all main.nf files in modules
NF_FILES=$(fd "^main.nf$" modules/)

template=".github/env-template.yml"

for file in $NF_FILES; do
    # Get the conda line
    conda_line=$(rg "^    conda" "$file")
    # Pull out the conda packages
    conda_packages=$(echo "$conda_line" | sed 's/^.*"\(.*\)".*$/\1/g')
    # Write the template to the environment.yml
    cat $template > $(dirname "$file")/environment.yml

    env_yaml=$(dirname "$file")/environment.yml
    # Put the dependencies in a environment.yaml using yq
    for package in $conda_packages; do
        yq -i '.dependencies += ["'"$package"'"]' $env_yaml
    done

    # Change the conda line to point to the new environment.yml
    sed -i "s|^    conda.*|    conda '$env_yaml'|g" "$file"
done
