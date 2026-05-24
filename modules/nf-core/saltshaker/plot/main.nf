process SALTSHAKER_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e93d703b195dd27cd920cee46669d3f51043216c12fd05168c937e93adf170e8/data':
        'community.wave.seqera.io/library/pip_saltshaker:e08e38a6d45f8f32' }"

    input:
    tuple val(meta), path(classify)

    output:
    tuple val(meta), path("*saltshaker.png"), emit: plot
    tuple val("${task.process}"), val('saltshaker'), val("1.0.0"), topic: versions, emit: versions_saltshaker

    when:
    task.ext.when == null || task.ext.when

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
