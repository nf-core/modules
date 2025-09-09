#!/usr/bin/env cached-nix-shell
#! nix-shell -p fd yq-go

# Find all the environment.ymls
ENV_FILES=$(fd "^environment.yml$" modules/)

for file in $ENV_FILES; do
    # Get the "name" field from the meta.yml next to the file
    name=$(yq -r '.name' $(dirname "$file")/meta.yml)
    # Add a name field in the enviroment.yml at the top of the file
    sed -i "1s|^|name: $name\n|" "$file"
done
