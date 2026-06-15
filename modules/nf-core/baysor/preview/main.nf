process BAYSOR_PREVIEW {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // images should be updated with wave containers when the issue of build number change is resolved
    container "quay.io/khersameesh24/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts), path(config)

    output:
    tuple val(meta), path("${prefix}_preview.html"), emit: html
    tuple val("${task.process}"), val('baysor'), eval("baysor --version"), topic: versions, emit: versions_baysor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export JULIA_NUM_THREADS=${task.cpus}

    baysor \\
        preview \\
        ${transcripts} \\
        --config ${config} \\
        --output ${prefix}_preview.html \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch "${prefix}_preview.html"
    """
}