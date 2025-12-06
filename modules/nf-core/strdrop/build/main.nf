process STRDROP_BUILD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4fc3114f6d67fadb826981d29ae67f8564ee81283184893c2db677d919b5b32d/data':
        'community.wave.seqera.io/library/pip_strdrop:df8d5dc993ea6848' }"

    input:
    tuple val(meta), path(training_set, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val('strdrop'), eval("strdrop --version | sed 's/.* //g'"), topic: versions, emit: versions_strdrop

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    strdrop \\
        build \\
        --training-set ./input \\
        ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    """
}
