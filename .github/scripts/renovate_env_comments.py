#!/usr/bin/env python3

import os
import re

# Define the directories to search for environment.yml files
directories = [
    "modules/nf-core/bowtie/align",
    "modules/nf-core/bowtie/build",
    "modules/nf-core/samtools/view",
    # Add more directories as needed
]

# Define a regex pattern to match dependencies
dependency_pattern = re.compile(r"^\s*-\s*(?P<channel>[\w-]+)::(?P<package>[\w-]+)=?(?P<version>[\d\.]*)")


# Function to add or verify renovate comments
def process_environment_file(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    updated_lines = []
    for line in lines:
        match = dependency_pattern.match(line)
        if match:
            channel = match.group("channel")
            package = match.group("package")
            version = match.group("version")
            comment = f"# renovate: datasource=conda depName={channel}/{package}\n"

            # Check if the comment is already present
            if comment not in updated_lines:
                updated_lines.append(comment)

        updated_lines.append(line)

    # Write the updated lines back to the file
    with open(file_path, "w") as file:
        file.writelines(updated_lines)


# Traverse directories and process each environment.yml file
for directory in directories:
    for root, _, files in os.walk(directory):
        for file in files:
            if file == "environment.yml":
                file_path = os.path.join(root, file)
                process_environment_file(file_path)
                print(f"Processed {file_path}")
