process SAMPLESHEETPARSER_CONVERT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:1.2.0--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheet)
    val target_version

    output:
    tuple val(meta), path("*.converted.csv"), emit: samplesheet
    tuple val("${task.process}"), val('samplesheet-parser'), eval("samplesheet --version | sed 's/samplesheet-parser //'"), topic: versions, emit: versions_samplesheet_parser

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!['v1', 'v2'].contains(target_version.toLowerCase())) {
        error "target_version must be 'v1' or 'v2', got: ${target_version}"
    }
    """
    samplesheet convert \\
        --to ${target_version} \\
        --output ${prefix}.converted.csv \\
        ${args} \\
        ${samplesheet}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.converted.csv
    """
}
