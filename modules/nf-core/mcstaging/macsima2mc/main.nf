// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MCSTAGING_MACSIMA2MC {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f571113c86e3216bf5fd714885f5521b813c4599a2852a15867d5bf3137a235/data' :
        'community.wave.seqera.io/library/python_macsima2mc:e3c1eb940a13b88b' }"

    input:
    tuple val(meta), path(input_dir), val(output_dir)

    output:
    tuple val(meta), path("${output_dir}/*")    , emit: out_dir
    tuple val("${task.process}"), val('macsima2mc'), eval('python -m pip show macsima2mc | grep "Version" | sed -e "s/Version: //g"'), topic: versions, emit: versions_macsima2mc

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "macsima2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }

    def args   = task.ext.args   ?: ''

    """
    mkdir -p ${output_dir}

    # iterate over individual cycle folders and create a single tiff file for each cycle
    for cycle in ${input_dir}/*Cycle*;
    do
        cycle_folder=\$(basename "\$cycle")
        macsima2mc \\
            -i \$cycle \\
            -o ${output_dir} \\
            ${args}
    done
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "macsima2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }

    """
    mkdir ${output_dir}
    mkdir ${output_dir}/well-rack-roi-exp
    touch ${output_dir}/well-rack-roi-exp/markers.csv
    mkdir ${output_dir}/well-rack-roi-exp/raw
    touch ${output_dir}/well-rack-roi-exp/raw/well-rack-roi-exp.ome.tif
    """
}
