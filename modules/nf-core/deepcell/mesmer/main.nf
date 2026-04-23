process DEEPCELL_MESMER {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/deepcell_mesmer:0.4.1_noentry"

    input:
    tuple val(meta) , path(img)
    tuple val(meta2), path(membrane_img)

    output:
    tuple val(meta), path("*.tif"), emit: mask
    tuple val("${task.process}"), val('deepcell'), eval("pip show 'DeepCell' | sed '2!d;s/Version: //'"), emit: versions_deepcell, topic: versions
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPCELL_MESMER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def membrane_command = membrane_img ? "--membrane-image $membrane_img" : ""

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "DEEPCELL_MESMER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tif
    """
}
