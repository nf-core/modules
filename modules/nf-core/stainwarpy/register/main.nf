process STAINWARPY_REGISTER {
    tag "$meta.id"
    label 'process_single'

    container "community.wave.seqera.io/library/pip_stainwarpy:d2ff09a87c79896b"

    input:
    tuple val(meta), path(hne_img)
    tuple val(meta2), path(multiplx_img)
    val fixed_img
    val final_sz

    output:
    tuple val(meta), path("0_final_channel_image.ome.tif")             , emit: reg_image
    tuple val(meta), path("registration_metrics_tform_map.json")       , emit: reg_metrics_tform
    tuple val(meta), path("feature_based_transformation_map.npy")      , emit: tform_map
    tuple val("${task.process}"), val('stainwarpy'), eval("stainwarpy --version | sed 's/.* //'"), emit: versions_stainwarpy_register, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    stainwarpy \\
        register \\
        ${multiplx_img} \\
        ${hne_img} \\
        . \\
        ${fixed_img} \\
        ${final_sz} \\
        ${args}
    """

    stub:

    """
    touch 0_final_channel_image.ome.tif
    touch registration_metrics_tform_map.json
    touch feature_based_transformation_map.npy
    """
}
