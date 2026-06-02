process BAYSOR_SEGFREE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // images should be updated with wave containers when the issue of build number change is resolved
    container "quay.io/khersameesh24/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts), path(config)

    output:
    tuple val(meta), path("${prefix}/ncvs.loom"), emit: ncvs
    tuple val("${task.process}"), val('baysor'), eval("baysor --version"), topic: versions, emit: versions_baysor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export JULIA_NUM_THREADS=${task.cpus}

    mkdir -p ${prefix}

    baysor \\
        segfree \\
        ${transcripts} \\
        --config ${config} \\
        --output ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/ncvs.loom"
    """
}