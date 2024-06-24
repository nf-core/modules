process DEEPCELL_MESMER {
    tag "$meta.id"
    label 'process_low'

    container "nf-core/deepcell_mesmer:0.4.1_noentry"

    input:
    tuple val(meta) , path(img)
    tuple val(meta2), path(membrane_img)

    // Output a .tif image, don't touch versions
    output:
    tuple val(meta), path("*.tif"), emit: mask
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def membrane_command = membrane_img ? "--membrane-image $membrane_img" : ""
    def VERSION = "0.4.1"

    """
    python /usr/src/app/run_app.py mesmer \\
        --squeeze \\
        --nuclear-image $img \\
        --output-directory . \\
        --output-name ${prefix}.tif \\
        $membrane_command \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepcell_mesmer: $VERSION
    END_VERSIONS
    """
}
