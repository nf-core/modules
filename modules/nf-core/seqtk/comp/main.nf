
process SEQTK_COMP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--h577a1d6_3':
        'biocontainers/seqtk:1.4--h577a1d6_3' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.seqtk_stats.tsv"), emit: seqtk_stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqtk comp \\
        ${args} \\
        ${fastx} > ${prefix}.seqtk_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk |& sed '/Version/!d; s/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" > ${prefix}.seqtk_stats.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk |& sed '/Version/!d; s/.* //')
    END_VERSIONS
    """
}
