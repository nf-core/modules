process VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/nf-core/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(segmentation_params), val(tile_index)
    path(algorithm_json) // Not passed as an arg; defined in segmentation parameters
    path(custom_weights) // Optional; also defined in segmentation parameters

    output:
    tuple val(meta), path("${prefix}/result_tiles/*.parquet")                                                                 , emit: segmented_tile
    tuple val("${task.process}"), val('vpt'), eval("pip show vpt | sed -n 's/Version: //p'")                                 , emit: versions_vpt, topic: versions
    tuple val("${task.process}"), val('vpt-plugin-cellpose2'), eval("pip show vpt-plugin-cellpose2 | sed -n 's/Version: //p'"), emit: versions_vptplugincellpose2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    vpt --verbose \\
        run-segmentation-on-tile \\
        $args \\
        --input-segmentation-parameters $segmentation_params \\
        --tile-index $tile_index
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/result_tiles
    touch ${prefix}/result_tiles/cell_0.parquet
    """
}
