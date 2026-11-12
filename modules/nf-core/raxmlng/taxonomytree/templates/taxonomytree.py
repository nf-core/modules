#!/usr/bin/env python3
import sys

prefix = "${meta.id}" if "${task.ext.prefix}" == "null" else "${task.ext.prefix}"


def parse_taxonomy(filepath):
    tab = chr(9)
    entries = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            seq_name, _, tax_str = line.partition(tab)
            ranks = [r.strip() for r in tax_str.split(";") if r.strip()]
            entries.append((seq_name, ranks))
    return entries


def build_trie(entries):
    root = {"children": {}, "leaves": []}
    for seq_name, ranks in entries:
        node = root
        for rank in ranks:
            if rank not in node["children"]:
                node["children"][rank] = {"children": {}, "leaves": []}
            node = node["children"][rank]
        node["leaves"].append(seq_name)
    return root


def node_to_newick(node):
    parts = []
    for child in node["children"].values():
        parts.append(node_to_newick(child))
    parts.extend(node["leaves"])
    if len(parts) == 1:
        return parts[0]
    return "(" + ",".join(parts) + ")"


def main(taxonomy_file, output_file):
    entries = parse_taxonomy(taxonomy_file)
    if not entries:
        sys.exit("No entries found in: " + taxonomy_file)
    root = build_trie(entries)
    newick = node_to_newick(root) + ";"
    with open(output_file, "w") as fh:
        print(newick, file=fh)


main("${taxonomy}", prefix + ".guide.nwk")

with open("versions.yml", "w") as fh:
    print('"${task.process}":', file=fh)
    print("    python: " + sys.version.split()[0], file=fh)
