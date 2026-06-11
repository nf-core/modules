process FINALETOOLKIT_AGGBW {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34f01d128ed135aedc33ddb62fced3911bef6d1a909694291b7184bf83719402/data'
        : 'community.wave.seqera.io/library/finaletoolkit:0.11.1--8fe5ba6ec9e2ec95'}"

    input:
    tuple val(meta), path(bigwig)
    tuple val(meta2), path(intervals)

    output:
    tuple val(meta), path("${prefix}.wig"), emit: wig
    tuple val("${task.process}"), val('finaletoolkit'), eval("finaletoolkit --version | sed 's/FinaleToolkit //g'"), topic: versions, emit: versions_finaletoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_agg"
    """
    finaletoolkit \\
        agg-bw \\
        ${bigwig} \\
        ${intervals} \\
        ${args} \\
        -o "${prefix}.wig"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_agg"
    """
    echo ${args}

    touch "${prefix}.wig"
    """
}
