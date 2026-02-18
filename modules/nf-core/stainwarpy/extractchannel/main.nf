process STAINWARPY_EXTRACTCHANNEL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c020e2c9696127244a6a2bbb1d945b9d78c18d3595bac54838e26bddcc87f521/data' :
        'community.wave.seqera.io/library/stainwarpy:0.2.4--c8bf19657f01e47a'}"

    input:
    tuple val(meta), path(multiplx_img)

    output:
    tuple val(meta), path("*_multiplexed_single_channel_img.ome.tif")                            , emit: single_ch_image
    tuple val("${task.process}"), val('stainwarpy'), eval("stainwarpy --version | sed 's/.* //'"), emit: versions_stainwarpy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stainwarpy \\
        extract-channel \\
        ${multiplx_img} \\
        . \\
        ${args}

    mv multiplexed_single_channel_img.ome.tif ${prefix}_multiplexed_single_channel_img.ome.tif
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_multiplexed_single_channel_img.ome.tif
    """
}
