process VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/nf-core/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(segmentation_params)
    path(algorithm_json) // Not passed as an arg but accessed by tool
    path(segmentation_tiles)

    output:
    tuple val(meta), path("${prefix}/*_mosaic_space.parquet")                                                                                , emit: mosaic_space
    tuple val(meta), path("${prefix}/*_micron_space.parquet")                                                                                , emit: micron_space
    tuple val("${task.process}"), val('vpt'), eval("pip show vpt | sed -n 's/Version: //p'")                                                 , emit: versions_vpt, topic: versions
    tuple val("${task.process}"), val('vpt-plugin-cellpose2'), eval("pip show vpt-plugin-cellpose2 | sed -n 's/Version: //p'")                , emit: versions_vptplugincellpose2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/result_tiles
    for segment in ${segmentation_tiles}; do
        cp -d \$segment ${prefix}/result_tiles
    done

    vpt --verbose \\
        compile-tile-segmentation \\
        ${args} \\
        --input-segmentation-parameters ${segmentation_params}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/segmented_micron_space.parquet
    touch ${prefix}/segmented_mosaic_space.parquet
    """

}
