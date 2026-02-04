process HMMCOPY_GENERATEMAP {
    tag "$fasta"
    label 'process_long'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmcopy:0.1.1--h2e03b76_7':
        'biocontainers/hmmcopy:0.1.1--h2e03b76_7' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    tuple val("${task.process}"), val('hmmcopy'), eval("echo 0.1.1"), topic: versions, emit: versions_hmmcopy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # build required indexes
    generateMap.pl -b \
        $args \
        $fasta

    # run
    generateMap.pl \
        $args \
        $fasta
    """
    stub:
    """
    touch ${fasta}.map.bw
    """
}
