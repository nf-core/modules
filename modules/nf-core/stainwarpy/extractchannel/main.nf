process STAINWARPY_EXTRACTCHANNEL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/stainwarpy:0.2.3--5966e23f2f7d254a' :
        'community.wave.seqera.io/library/stainwarpy:0.2.3--2c8b18a5e6d93e4a'}"

    input:
    tuple val(meta), path(multiplx_img)

    output:
    tuple val(meta), path("multiplexed_single_channel_img.ome.tif")    , emit: single_ch_image
    tuple val("${task.process}"), val('stainwarpy'), eval("stainwarpy --version | sed 's/.* //'"), emit: versions_stainwarpy_extractchannel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    stainwarpy \\
        extract-channel \\
        ${multiplx_img} \\
        . \\
        ${args}
    """

    stub:

    """
    touch multiplexed_single_channel_img.ome.tif
    """
}
