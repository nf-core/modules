process UCSC_BEDTOBIGBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bedtobigbed:447--h954228d_0' :
        'biocontainers/ucsc-bedtobigbed:447--h954228d_0' }"

    input:
    tuple val(meta), path(bed)
    path  sizes
    path  autosql

    output:
    tuple val(meta), path("*.bigBed"), emit: bigbed
    tuple val("${task.process}"), val('ucsc'), val('477'), topic: versions, emit: versions_ucsc
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
