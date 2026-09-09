process BAYSOR_PREVIEW {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f716efdfa817ee57bf59b78248f672b503807bfb8608ad3d3976f2e706ca9fb4/data' :
        'community.wave.seqera.io/library/baysor:0.7.1--6fd896e03359bae6'}"

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
