process SAMPLESHEETPARSER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:1.2.0--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheets, stageAs: "?/*")
    val target_version

    output:
    tuple val(meta), path("*.merged.csv"), emit: samplesheet
    tuple val(meta), path("*.merge.json"), emit: json
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
    samplesheet merge \\
        --to ${target_version} \\
        --output ${prefix}.merged.csv \\
        --format json \\
        --force \\
        ${args} \\
        ${samplesheets} > ${prefix}.merge.json || true
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.csv
    echo '{"has_conflicts": false, "sample_count": 0, "conflicts": [], "warnings": []}' > ${prefix}.merge.json
    """
}
