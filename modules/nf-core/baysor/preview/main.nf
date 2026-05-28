process BAYSOR_PREVIEW {
    tag "${meta.id}"
    label 'process_medium'

    container "khersameesh24/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts), path(config)

    output:
    tuple val(meta), path("${prefix}/preview.html"), emit: preview_html
    tuple val("${task.process}"), val('baysor'), eval("baysor --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo unknown"), topic: versions, emit: versions_baysor

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BAYSOR_PREVIEW module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export JULIA_NUM_THREADS=${task.cpus}

    mkdir -p ${prefix}

    baysor preview \\
    ${transcripts} \\
    --config ${config} \\
    --output ${prefix}/preview.html
    ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BAYSOR_PREVIEW module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/preview.html
    """
}