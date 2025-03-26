process MCSTAGING_MACSIMA2MC {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/macsima2mc:v1.2.13"

    input:
    tuple val(meta), path(input_dir), val(output_dir)

    output:
    tuple val(meta), path("${output_dir}/*")    , emit: out_dir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "macsima2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def VERSION = '1.2.13'

    """
    python /staging/macsima2mc/macsima2mc.py \
        -i ${input_dir} \
        -o ${output_dir} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsima2mc: ${VERSION}
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "macsima2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }

    def VERSION = '1.2.13'

    """
    mkdir input_dir
    mkdir output_dir
    mkdir output_dir/well-rack-roi-exp
    touch output_dir/markers.csv
    mkdir output_dir/well-rack-roi-exp/raw
    touch output_dir/well-rack-roi-exp/raw/well-rack-roi-exp.ome.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macsima2mc: ${VERSION}
    END_VERSIONS
    """
}
