process NUMORPH_RESAMPLE {
    tag "$meta.id"
    label 'process_low'

    container "nf-core/numorph_analyze:1.0.0"

    input:
    tuple val(meta), path(stitch_directory), path(parameter_file)

    output:
    tuple val(meta), path("results/resampled/")                        , emit: resampled
    tuple val("${task.process}"), val('numorph_resample'), val('1.0.0'), emit: versions_numorph_analyze, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    mkdir -p ./results/stitched

    ln -sr ${stitch_directory}/* results/stitched

    # resolve symlinks and paths
    stitch_directory=\$(readlink -f ./results/stitched)
    parameter_file=\$(readlink -f ${parameter_file})
    results_dir=\$(readlink -f ./results)

    numorph_analyze 'input_dir' \$stitch_directory \\
        'output_dir' \$results_dir \\
        'parameter_file' \$parameter_file \\
        'sample_name' ${meta.id} \\
        'stage' 'resample' \\
        'NM_variables' '' \\
        'use_processed_images' 'stitched'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p results/resampled

    touch results/resampled/${prefix}_resampled.tif
    """
}
