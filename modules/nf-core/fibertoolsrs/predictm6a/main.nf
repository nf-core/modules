process FIBERTOOLSRS_PREDICTM6A {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fibertools-rs:0.7.1--h3b373d1_0':
        'biocontainers/fibertools-rs:0.7.1--h3b373d1_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_m6a"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ft \\
        predict-m6a \\
        $args \\
        --threads ${task.cpus} \\
        $bam \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fibertools-rs: \$(ft --version | sed -E 's/.* ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_m6a"

    """
    echo $args

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fibertools-rs: \$(ft --version | sed -E 's/.* ([0-9]+\\.[0-9]+\\.[0-9]+).*/\\1/' )
    END_VERSIONS
    """
}
