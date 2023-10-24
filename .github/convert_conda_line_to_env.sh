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
    # Put the dependencies in a environment.yaml using yq
    yq -i '.dependencies += ["'"$conda_packages"'"]' $(dirname "$file")/environment.yml
done
