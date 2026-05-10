process SAMPLESHEETPARSER_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:1.2.0--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheet)

    output:
    tuple val(meta), path("*.filtered.csv"), emit: samplesheet, optional: true
    tuple val(meta), path("*.filter.json"),  emit: json
    tuple val("${task.process}"), val('samplesheet-parser'), eval("samplesheet --version | sed 's/samplesheet-parser //'"), topic: versions, emit: versions_samplesheet_parser

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samplesheet filter \\
        --output ${prefix}.filtered.csv \\
        --format json \\
        ${args} \\
        ${samplesheet} > ${prefix}.filter.json || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.csv
    echo '{"matched_count": 0, "total_count": 0, "output_path": null, "criteria": {}}' > ${prefix}.filter.json
    """
}
