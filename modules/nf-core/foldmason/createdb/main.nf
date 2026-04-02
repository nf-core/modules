process FOLDMASON_CREATEDB {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc0b27a73277c851f0a1df09fc5fa4f0ff80832e6779f143420a7e2321ac62b4/data' :
            'community.wave.seqera.io/library/foldmason:4.dd3c235--375019dcc3c5fe6f' }"

    input:
    tuple val(meta) , path(structures)

    output:
    tuple val(meta), path("${prefix}*"), emit: db
    tuple val("${task.process}"), val('foldmason'), eval('foldmason version'), emit: versions_foldmason, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    foldmason createdb \\
        ${structures} \\
        ${prefix} \\
        $args \\
        --threads $task.cpus
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    """
}
