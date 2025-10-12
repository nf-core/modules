process VIZGENPOSTPROCESSING_COMPILETILESEGMENTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'nf-core/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(segmentation_params)
    path(algorithm_json) // Not passed as an arg but accessed by tool
    path(segmentation_tiles)

    output:
    tuple val(meta), path("${prefix}/*_mosaic_space.parquet"), emit: mosaic_space
    tuple val(meta), path("${prefix}/*_micron_space.parquet"), emit: micron_space
    path  "versions.yml"                                     , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/segmented_micron_space.parquet
    touch ${prefix}/segmented_mosaic_space.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """

}
