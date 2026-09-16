process CUSTOM_TAXONOMYTREE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // quay.io/biocontainers/python:3.11 has no build-hash-suffixed tag to pin to (unlike
    // real bioconda-recipe images); pin by digest instead so the underlying image can't
    // silently drift between runs.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python@sha256:b322907f8e52b2055ccad4e46848d28a4a5631b403116cc80ddf61ec8601e05e' }"

    input:
    tuple val(meta), path(taxonomy)

    output:
    tuple val(meta), path("*.guide.nwk"), emit: guide_tree
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 - "${taxonomy}" "${prefix}.guide.nwk" << 'PYEOF'
import sys

def parse_taxonomy(filepath):
    tab = chr(9)
    entries = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            seq_name, _, tax_str = line.partition(tab)
            ranks = [r.strip() for r in tax_str.split(';') if r.strip()]
            entries.append((seq_name, ranks))
    return entries

def build_trie(entries):
    root = {'children': {}, 'leaves': []}
    for seq_name, ranks in entries:
        node = root
        for rank in ranks:
            if rank not in node['children']:
                node['children'][rank] = {'children': {}, 'leaves': []}
            node = node['children'][rank]
        node['leaves'].append(seq_name)
    return root

def node_to_newick(node):
    parts = []
    for child in node['children'].values():
        parts.append(node_to_newick(child))
    parts.extend(node['leaves'])
    if len(parts) == 1:
        return parts[0]
    return '(' + ','.join(parts) + ')'

def main(taxonomy_file, output_file):
    entries = parse_taxonomy(taxonomy_file)
    if not entries:
        sys.exit('No entries found in: ' + taxonomy_file)
    root = build_trie(entries)
    newick = node_to_newick(root) + ';'
    with open(output_file, 'w') as fh:
        print(newick, file=fh)

main(sys.argv[1], sys.argv[2])
PYEOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.guide.nwk
    """
}
