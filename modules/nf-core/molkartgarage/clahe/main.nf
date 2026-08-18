process MOLKARTGARAGE_CLAHE {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/schapirolabor/molkart-local:v0.1.1"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("*.tiff") , emit: img_clahe
    tuple val("${task.process}"), val('clahe'), eval("python /local/scripts/molkart_clahe.py --version"), emit: versions_clahe, topic: versions
    tuple val("${task.process}"), val('scikit-image'), eval("python -c 'import skimage; print(skimage.__version__)'"), emit: versions_scikitimage, topic: versions

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

    sed -i -E 's/UUID="urn:uuid:[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}"/                                                    /g' ${prefix}.tiff
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tiff
    """
}
