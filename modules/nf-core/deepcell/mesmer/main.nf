process DEEPCELL_MESMER {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/deepcell_mesmer:0.4.1_noentry"

    input:
    tuple val(meta) , path(img)
    tuple val(meta2), path(membrane_img)

    output:
    tuple val(meta), path("*.tif"), emit: mask
    tuple val("${task.process}"), val('deepcell_mesmer'), val("0.4.1"), emit: versions_mesmer, topic: versions
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def membrane_command = membrane_img ? "--membrane-image $membrane_img" : ""
    def VERSION          = "0.4.1"

    """
    python /usr/src/app/run_app.py mesmer \\
        --squeeze \\
        --nuclear-image $img \\
        --output-directory . \\
        --output-name ${prefix}.tif \\
        $membrane_command \\
        $args
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tif
    """
}
