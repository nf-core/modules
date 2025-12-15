process STRDROP_BUILD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f860e6cdc0d4222f89145d5e5f6aba15368eefc50b65bc78890613d976344a7f/data':
        'community.wave.seqera.io/library/pip_strdrop:b1aa6c1a4a3357f2' }"

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
