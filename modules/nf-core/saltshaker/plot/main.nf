process SALTSHAKER_PLOT {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5a902cc9f161d602fde9c268a509be2f593cfac7ed4cdc2219f630e02e43b2ec/data':
        'community.wave.seqera.io/library/pip_saltshaker:be40ca61bbf77cf2' }"

    input:
    tuple val(meta), path(classify)

    output:
    tuple val(meta), path("*saltshaker.png"), emit: plot
    tuple val("${task.process}"), val('saltshaker'), val("1.1.1"), topic: versions, emit: versions_saltshaker


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    saltshaker plot \\
        --prefix $prefix \\
        --input-dir . \\
        $args

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.saltshaker.png
    """
}
