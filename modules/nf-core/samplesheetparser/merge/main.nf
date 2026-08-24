process SAMPLESHEETPARSER_MERGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samplesheet-parser:2.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/samplesheet-parser:2.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(samplesheets)
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
    if (!['v1', 'v2'].contains(target_version.toLowerCase())) {
        error "target_version must be 'v1' or 'v2', got: ${target_version}"
    }
    """
    samplesheet merge \\
        --to ${target_version} \\
        --output ${prefix}.merged.csv \\
        --format json \\
        ${args} \\
        ${samplesheets} > ${prefix}.merge.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.csv
    touch ${prefix}.merge.json
    """
}
