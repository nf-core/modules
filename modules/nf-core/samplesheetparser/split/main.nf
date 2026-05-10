process SAMPLESHEETPARSER_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:1.2.0--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheet)
    val by

    output:
    tuple val(meta), path("split/*.csv"), emit: samplesheets
    tuple val(meta), path("*.split.json"), emit: json
    tuple val("${task.process}"), val('samplesheet-parser'), eval("samplesheet --version | sed 's/samplesheet-parser //'"), topic: versions, emit: versions_samplesheet_parser

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!['project', 'lane'].contains(by.toLowerCase())) {
        error "by must be 'project' or 'lane', got: ${by}"
    }
    """
    mkdir -p split

    samplesheet split \\
        --by ${by} \\
        --output-dir split \\
        --format json \\
        ${args} \\
        ${samplesheet} > ${prefix}.split.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p split
    touch split/stub_SampleSheet.csv
    echo '{"by": "${by}", "files": {}, "sample_counts": {}, "warnings": []}' > ${prefix}.split.json
    """
}
