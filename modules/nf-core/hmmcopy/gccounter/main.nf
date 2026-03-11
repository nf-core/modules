process HMMCOPY_GCCOUNTER {
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7' :
        'biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.wig"), emit: wig
    tuple val("${task.process}"), val('hmmcopy'), eval("echo 0.1.1"), topic: versions, emit: versions_hmmcopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_gc"
    """
    gcCounter \
        $args \
        ${fasta} > ${prefix}.wig
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_gc"
    """
    touch ${prefix}.wig
    """
}
