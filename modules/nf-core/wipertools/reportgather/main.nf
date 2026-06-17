process WIPERTOOLS_REPORTGATHER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wipertools:1.1.5--pyhdfd78af_0':
        'quay.io/biocontainers/wipertools:1.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("${prefix}.report"), emit: gathered_report
    tuple val("${task.process}"), val('wipertools'), eval("wipertools reportgather --version"), topic: versions, emit: versions_wipertools


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (report.any { file -> file.name == "${prefix}.report" }) {
        error 'Output file name "${prefix}.report" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    wipertools \\
        reportgather \\
        -r ${report} \\
        -f ${prefix}.report \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_gather"

    // Check if the output file name is in the list of input files
    if (report.any { file -> file.name == "${prefix}.report" }) {
        error 'Output file name "${prefix}.report" matches one of the input files. Use \"task.ext.prefix\" to disambiguate!.'
    }

    """
    touch ${prefix}.report
    """
}
