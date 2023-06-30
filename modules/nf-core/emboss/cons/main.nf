process EMBOSS_CONS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::emboss=6.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--h86d058a_5':
        'biocontainers/emboss:6.6.0--h86d058a_5' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.")   , emit: consensu
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = (args =~ "name") ? '' : "-name $prefix"
    """
    cons \\
        ${args} \\
        ${name} \\
        -sequence $fasta \\
        -outseq ${prefix} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(cons --version 2>&1) | sed 'EMBOSS://' ))
    END_VERSIONS
    """
}
