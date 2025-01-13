process WIPERTOOLS_REPORTGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.3--pyhdfd78af_0':
        'biocontainers/wipertools:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(report)   // channel: [ val(meta), [ *.report ] ]

    output:
    tuple val(meta), path("${prefix}.report"), emit: gathered_report    // channel: [ val(meta), *_gather.report ]
    path "versions.yml"                      , emit: versions           // channel: [ versions.yml ]

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (report.any { it.name == "${prefix}.report" }) {
        error 'Output file name "${prefix}.report}" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    wipertools \\
        reportgather \\
        -r $report \\
        -f ${prefix}.report \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (report.any { it.name == "${prefix}.report" }) {
        error 'Output file name "${prefix}.report}" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    touch ${prefix}.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wipertools reportgather: \$(wipertools reportgather --version)
    END_VERSIONS
    """
}
