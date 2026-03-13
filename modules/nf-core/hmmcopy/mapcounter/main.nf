process HMMCOPY_MAPCOUNTER {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7'
        : 'biocontainers/hmmcopy:0.1.1--h2e03b76_7'}"

    input:
    tuple val(meta), path(bigwig)

    output:
    tuple val(meta), path("*.wig"), emit: wig
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('hmmcopy'), val("0.1.1"), topic: versions, emit: versions_hmmcopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_map"
    """
    mapCounter \
        ${args} \
        ${bigwig} > ${prefix}.wig
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_map"
    """
    touch ${prefix}.wig
    """
}
