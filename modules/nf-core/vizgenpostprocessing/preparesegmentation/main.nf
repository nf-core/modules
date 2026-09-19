process VIZGENPOSTPROCESSING_PREPARESEGMENTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/nf-core/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(um_to_mosaic_file)
    path(algorithm_json)
    val(images_regex)

    output:
    tuple val(meta), path("${prefix}/*.json")                                                                                 , emit: segmentation_files
    tuple val("${task.process}"), val('vpt'), eval("pip show vpt | sed -n 's/Version: //p'")                                 , emit: versions_vpt, topic: versions
    tuple val("${task.process}"), val('vpt-plugin-cellpose2'), eval("pip show vpt-plugin-cellpose2 | sed -n 's/Version: //p'"), emit: versions_vptplugincellpose2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    vpt --verbose \\
        prepare-segmentation \\
        ${args} \\
        --segmentation-algorithm ${algorithm_json} \\
        --input-images "${input_images}/${images_regex}" \\
        --input-micron-to-mosaic ${um_to_mosaic_file} \\
        --output-path ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/algorithm_specification.json
    """
}
