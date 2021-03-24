#!/usr/bin/env python
import os
from textwrap import dedent

results = {}
version_files = [x for x in os.listdir(".") if x.endswith(".version.txt")]
for version_file in version_files:

    software = version_file.replace(".version.txt", "")

    with open(version_file) as fin:
        version = fin.read().strip()

    results[software] = version

# Dump to YAML
print(
    dedent(
        """\
        id: 'software_versions'
        section_name: 'nf-core/rnaseq Software Versions'
        section_href: 'https://github.com/nf-core/rnaseq'
        plot_type: 'html'
        description: 'are collected at run time from the software output.'
        data: |
            <dl class="dl-horizontal">
        """
    )
)
for k, v in sorted(results.items()):
    print(f"        <dt>{k}</dt><dd><samp>{v}</samp></dd>")
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.tsv", "w") as f:
    for k, v in sorted(results.items()):
        f.write(f"{k}\t{v}\n")
