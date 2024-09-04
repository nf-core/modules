process MCSTAGING_IMC2MC {
    tag "$meta.id"
    label 'process_single'
    container "ghcr.io/schapirolabor/imc2mc:0.0.2"

    input:
    tuple val(meta) , path(txtfile)

    output:
    tuple val(meta), path("*.tif"), emit: tif
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "imc2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python /imc2mc/scripts/imc2mc.py \
        -i ${txtfile} \
        -o "${prefix}.tif" \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imc2mc: \$(python /imc2mc/scripts/imc2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "imc2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        imc2mc: \$(python /imc2mc/scripts/imc2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
