process EMBOSS_REVSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5':
        'biocontainers/emboss:6.6.0--h86d058a_5' }"

    input:
    tuple val(meta), path(sequences)

    output:
    tuple val(meta), path("*.${sequences.name - ~/.*\./}"), emit: revseq
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = sequences.name - ~/.*\./
    def outfile = "${prefix}.rev.${suffix}"
    """
    revseq \\
        $args \\
        $sequences \\
        $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(echo \$(revseq -version 2>&1) | sed 's/EMBOSS://')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = sequences.name - ~/.*\./
    def outfile = "${prefix}.rev.${suffix}"
    """
    touch ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(echo \$(revseq -version 2>&1) | sed 's/EMBOSS://')
    END_VERSIONS
    """
}
