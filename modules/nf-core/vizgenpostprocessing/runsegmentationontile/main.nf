process VIZGENPOSTPROCESSING_RUNSEGMENTATIONONTILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(segmentation_params), val(tile_index)
    path(algorithm_json) // Not passed as an arg; defined in segmentation parameters
    path(custom_weights) // Optional; also defined in segmentation parameters

    output:
    tuple val(meta), path("${prefix}/result_tiles/*.parquet"), emit: segmented_tile
    path  "versions.yml"                                     , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/result_tiles
    touch ${prefix}/result_tiles/cell_0.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """
}
