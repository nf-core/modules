#!/usr/bin/env sh

# Get all of the meta.yaml files
METAS=$(fd meta.yml)

# Define the output file path
output_file=".github/CODEOWNERS-tmp"

# Use yq to extract the "authors" array and convert it to a .gitignore format
for file in $METAS; do
    # Print the directory the file is in first
    path=$(echo "$file" | sed 's/\/meta.yml//')
    # Add a double star to the end of the path
    path="$path/**"

    authors=$(yq '.authors | .[]' "$file" | sed 's/^//')
    # Remove quotes from authors
    authors=$(echo "$authors" | sed 's/"//g')
    echo "$path" $authors >> $output_file
done

cat $output_file >> ".github/CODEOWNERS"
