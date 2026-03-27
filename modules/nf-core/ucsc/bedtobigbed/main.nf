process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:482--hdc0a859_0' :
        'biocontainers/ucsc-bedtobigbed:482--hdc0a859_0' }"

    input:
    tuple val(meta), path(bed)
    path  sizes
    path  autosql

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def as_option = autosql ? "-as=${autosql}" : ""
    """
    bedToBigBed \\
        $bed \\
        $sizes \\
        $as_option \\
        $args \\
        ${prefix}.bigBed

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bigBed
    """
}
