process STAINWARPY_TRANSFORMSEGMASK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c020e2c9696127244a6a2bbb1d945b9d78c18d3595bac54838e26bddcc87f521/data' :
        'community.wave.seqera.io/library/stainwarpy:0.2.4--c8bf19657f01e47a'}"

    input:
    tuple val(meta), path(hne_img)
    tuple val(meta2), path(multiplx_img)
    tuple val(meta3), path(seg_mask)
    tuple val(meta4), path(tform_map)
    val fixed_img
    val final_sz

    output:
    tuple val(meta), path("*_transformed_segmentation_mask.ome.tif")                             , emit: transformed_seg_mask
    tuple val("${task.process}"), val('stainwarpy'), eval("stainwarpy --version | sed 's/.* //'"), emit: versions_stainwarpy, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    stainwarpy \\
        transform-seg-mask \\
        ${seg_mask} \\
        ${multiplx_img} \\
        ${hne_img} \\
        . \\
        ${tform_map} \\
        ${fixed_img} \\
        ${final_sz} \\
        ${args}

    mv transformed_segmentation_mask.ome.tif ${prefix}_transformed_segmentation_mask.ome.tif
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_transformed_segmentation_mask.ome.tif
    """
}
