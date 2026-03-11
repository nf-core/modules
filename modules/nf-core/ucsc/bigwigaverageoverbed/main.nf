process UCSC_BIGWIGAVERAGEOVERBED {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bigwigaverageoverbed:482--h0b57e2e_0' :
        'biocontainers/ucsc-bigwigaverageoverbed:482--h0b57e2e_0' }"

    input:
    tuple val(meta), path(bed)
    path bigwig

    output:
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // BUG: bigWigAverageOverBed cannot handle ensembl seqlevels style
    """
    bigWigAverageOverBed \\
        $args \\
        $bigwig \\
        $bed \\
        ${prefix}.tab
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tab

    """
}
