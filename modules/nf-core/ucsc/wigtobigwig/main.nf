process UCSC_WIGTOBIGWIG {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/0394161be68e8dd5b30a47f0b19ffa00cb3226bb2e6c9fe3ec89e571a50b871d/data' :
        'community.wave.seqera.io/library/ucsc-wigtobigwig:482--7b910cc21c32327e' }"

    input:
    tuple val(meta), path(wig)
    path sizes

    output:
    tuple val(meta), path("${prefix}.bw"), emit: bw
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    wigToBigWig \\
        $args \\
        $wig \\
        $sizes \\
        ${prefix}.bw
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw
    """
}
