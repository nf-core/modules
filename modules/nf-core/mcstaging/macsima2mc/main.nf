process MCSTAGING_MACSIMA2MC {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/numpy_python_pip_macsima2mc:b42790f2c05a215a' :
        'community.wave.seqera.io/library/numpy_python_pip_macsima2mc:3fa6bf589964777e' }"

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
    macsima2mc \
        -i ${input_dir} \
        -o ${output_dir} \
        ${args}
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
