process WIPERTOOLS_REPORTGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("${prefix}.report"), emit: report_out
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    prefix   = prefix+"_gathered"
    """
    wipertools \\
        reportgather \\
        -r $reports \\
        -f ${prefix}.report \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    prefix   = prefix+"_gathered"
    """
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """
}
