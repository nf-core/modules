process MOLKARTGARAGE_CLAHE {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/schapirolabor/molkart-local:v0.1.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.tiff") , emit: img_clahe
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Molkartgarage/clahe module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /local/scripts/molkart_clahe.py \
        --input ${image} \
        --output ${prefix}.tiff \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        molkart_clahe: \$(python /local/scripts/molkart_clahe.py --version)
        scikit-image: 0.19.2
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        molkart_clahe: \$(python3 /local/scripts/molkart_clahe.py  --version)
        scikit-image: 0.19.2
    END_VERSIONS
    """
}
