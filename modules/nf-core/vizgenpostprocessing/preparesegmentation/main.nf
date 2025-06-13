process VIZGENPOSTPROCESSING_PREPARESEGMENTATION {
    tag "$meta.id"
    label 'process_single'

    container 'ghcr.io/wehi-soda-hub/vizgen-postprocessing_container:v0.1.1'

    input:
    tuple val(meta), path(input_images), path(um_to_mosaic_file)
    path(algorithm_json)
    val(images_regex)

    output:
    tuple val(meta), path("${prefix}/*.json"), emit: segmentation_files
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "VPT is unavailable via Conda. Please use Docker / Singularity / Apptainer / Podman instead."
    }
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "VPT is unavailable via Conda. Please use Docker / Singularity / Apptainer / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/algorithm_specification.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: \$( pip show vpt | grep Version | sed -e "s/Version: //g" )
        vpt-plugin-cellpose2: \$( pip show vpt-plugin-cellpose2 | grep Version | sed -e "s/Version: //g" )
    END_VERSIONS
    """
}
