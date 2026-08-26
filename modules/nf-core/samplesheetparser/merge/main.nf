process SAMPLESHEETPARSER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:2.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:2.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheets, stageAs: "input*/*")
    val target_version

    output:
    tuple val(meta), path("*.merged.csv"), emit: samplesheet
    tuple val(meta), path("*.merge.json"), emit: json
    tuple val("${task.process}"), val('samplesheet-parser'), eval("samplesheet --version | sed 's/samplesheet-parser //'"), topic: versions, emit: versions_samplesheetparser

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def to_norm = target_version.toLowerCase()
    def sheets  = samplesheets instanceof List ? samplesheets : [ samplesheets ]
    if (!['v1', 'v2'].contains(to_norm)) {
        error "target_version must be 'v1' or 'v2', got: ${target_version}"
    }
    if (sheets.size() < 2) {
        error "SAMPLESHEETPARSER_MERGE requires at least two input sheets, got ${sheets.size()}"
    }
    // `samplesheet merge` exits 1 when the report contains conflicts or warnings
    // (and 2 on hard errors). Tolerate exit 1 so a warning does not kill the task;
    // downstream code should inspect `has_conflicts` in the JSON report.
    """
    samplesheet merge \\
        --to ${to_norm} \\
        --output ${prefix}.merged.csv \\
        --format json \\
        ${args} \\
        ${samplesheets} > ${prefix}.merge.json || [ \$? -eq 1 ]
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.csv
    touch ${prefix}.merge.json
    """
}
