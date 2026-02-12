process RTGTOOLS_ROCPLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val(meta), path("*.svg"), emit: svg
    tuple val("${task.process}"), val('rgtools'), eval("rtg version | head -n 1 | sed 's/Product: RTG Tools //'"), topic: versions, emit: versions_rtgtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory.toGiga() + "G"
    """
    rtg RTG_MEM=${avail_mem} rocplot \\
        ${args} \\
        --png ${prefix}.png \\
        --svg ${prefix}.svg \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    touch ${prefix}.svg
    """
}
