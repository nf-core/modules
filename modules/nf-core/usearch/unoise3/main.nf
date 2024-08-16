process USEARCH_UNOISE3 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/usearch:12.0_beta--h9ee0642_1':
        'biocontainers/usearch:12.0_beta--h9ee0642_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    usearch \\
        -unoise3 $fasta \\
        $args \\
        -threads $task.cpus \\
        -zotus ${prefix}_features.fasta \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch| head -n 1 | cut -f 2 -d "v" |cut -f 1 -d " ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_features.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch| head -n 1 | cut -f 2 -d "v" |cut -f 1 -d " ")
    END_VERSIONS
    """
}
