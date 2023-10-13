process KALIGN_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kalign3=3.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kalign3:3.3.5--hdbdd923_0':
        'biocontainers/kalign3:3.3.5--hdbdd923_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln"), emit: alignment
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kalign \\
        $args \\
        -i $fasta \\
        -o ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kalign : \$(echo \$(kalign -v) | sed 's/kalign //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kalign : \$(echo \$(kalign -v) | sed 's/kalign //g' )
    END_VERSIONS
    """
}
