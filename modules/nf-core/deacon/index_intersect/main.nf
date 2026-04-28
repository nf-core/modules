process DEACON_INDEX_INTERSECT {
    tag "fasta"

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta), path(index)    // multiple deacon .idx index files to intersect

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: constructing a new meta.id would be more tricky here since there can be many input genomes
    // def meta_out = meta_index + [
    //     id: meta_index.id + '__' +
    //         meta_genome.collect{ it.id }.sort().join('_')
    // ]
    // def meta_out = meta + [
    //     id: meta.collect{ it.id }.sort().join('_')
    // ]

    """
    deacon \\
        index \\
        intersect \\
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
