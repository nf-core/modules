// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SAMPLESHEETPARSER_VALIDATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:2.6.0--pyhdfd78af_1' :
        'quay.io/biocontainers/samplesheet-parser:2.6.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(samplesheet)
    path(schema)

    output:
    tuple val(meta), path("*.validation.json"), emit: json
    tuple val("${task.process}"), val('samplesheet-parser'), eval("samplesheet --version | sed 's/samplesheet-parser //'"), topic: versions, emit: versions_samplesheetparser

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def schema_arg = schema ? "--schema ${schema}" : ''
    """
    samplesheet validate \\
        ${args} \\
        ${schema_arg} \\
        ${samplesheet} > ${prefix}.validation.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.validation.json
    """
}
