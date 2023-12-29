#!/usr/bin/env bash

# Get all of the meta.yaml files
METAS=$(fd meta.yml)

# Define the output file path
output_file=".github/CODEOWNERS-tmp"
rm $output_file

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

# Generate it from scratch
cat ".github/manual_CODEOWNERS" > ".github/CODEOWNERS"
# Remove duplicate lines and then sort and only print the line not the common part
cat "$output_file" | sort | uniq  >> ".github/CODEOWNERS"
