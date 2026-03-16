process NUMORPH_INTENSITY {
    tag "$meta.id"
    label 'process_low'

    container "nf-core/numorph_preprocessing:1.0.0"


    input:
    tuple val(meta), path(img_directory), path(parameter_file)

    output:
    tuple val(meta), path("results/variables/")                         , emit: variables
    tuple val(meta), path("results/samples/")                           , emit: samples
    tuple val(meta), path("results/NM_variables.mat")                   , emit: NM_variable

    tuple val("${task.process}"), val('numorph_intensity'), val('1.0.0'), emit: versions_numorph_intensity, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ./results

    # resolve symlinks and paths
    img_dir=\$(readlink -f ${img_directory})
    parameter_file=\$(readlink -f ${parameter_file})
    results_dir=\$(readlink -f ./results)


    numorph_preprocessing 'input_dir' \$img_dir 'output_dir' \$results_dir 'parameter_file' \$parameter_file 'sample_name' ${meta.id} 'stage' 'intensity'

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ./results/samples/intensity_adjustment
    mkdir -p ./results/variables

    touch ./results/variables/adj_params.mat
    touch ./results/variables/path_table.mat
    touch ./results/variables/thresholds.mat
    touch ./results/NM_variables.mat
    touch ./results/samples/intensity_adjustment/${meta.id}.png
    """
}
