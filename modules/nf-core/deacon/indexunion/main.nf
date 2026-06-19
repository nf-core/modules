process DEACON_INDEX_UNION {
    tag "fasta"

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta), path(index)    // multiple deacon .idx index files to combine

    // TODO: optionally could accept two inputs similar to index_diff:
    // the first would be the main index to start from, providing meta_index.id
    // the second would be a 1 or more indexes to combine with

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: constructing a new meta.id would be more tricky here since there can be many input genomes
    // def meta_out = meta + [
    //     id: meta.collect{ it.id }.sort().join('_')
    // ]

    """
    deacon \\
        index \\
        union \\
        ${args} \\
        ${index} \\
        > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.idx
    """
}
