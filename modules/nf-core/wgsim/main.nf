process WGSIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wgsim:1.0--h5bf99c6_4':
        'biocontainers/wgsim:1.0--h5bf99c6_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    path "versions.yml",              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    wgsim \\
        $args \\
        $fasta \\
        ${prefix}_R1.fastq \\
        ${prefix}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wgsim: \$(wgsim 2>&1 | grep "Version" | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_R1.fastq
    touch ${prefix}_R2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wgsim: \$(wgsim 2>&1 | grep "Version" | sed 's/Version: //')
    END_VERSIONS
    """
}
