process DEACON_INDEXDIFF {
    tag "fasta"

    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.13.2--h7ef3eeb_1':
        'biocontainers/deacon:0.13.2--h7ef3eeb_0' }"

    input:
    tuple val(meta_index), path(index)      // main deacon .idx index file
    tuple val(meta_genome), path(genome)    // a single fasta or .idx file to subtract from the main index

    output:
    tuple val(meta_index), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_index.id}.diff"

    """
    deacon \\
        index \\
        diff \\
        ${args} \\
        ${index} \\
        ${genome} \\
        > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta_index.id}"
    """
    touch ${prefix}.idx
    """
}
