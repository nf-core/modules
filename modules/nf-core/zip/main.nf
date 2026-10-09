process ZIP {
    tag "$prefix"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/171e3d1c631246d5afa0f01e78e0ff6a6ee78012dd17a7e8eccc871b65eac258/data' :
        'community.wave.seqera.io/library/p7zip:16.02--94efedc6cfea7db9' }"

    input:
    tuple val(meta), path(files, stageAs: "inputs/*")

    output:
    tuple val(meta), path("${prefix}.zip"), emit: zipped_archive
    tuple val("${task.process}"), val('zip'), eval("7za i 2>&1 | grep 'p7zip Version' | sed 's/.*p7zip Version //; s/(.*//'"), emit: versions_zip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : 'zipped_files')
    """
    7z \\
        a \\
        -l \\
        $args \\
        "${prefix}.zip" ./inputs/*
    """

    stub:
    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : 'zipped_files')
    """
    touch ${prefix}.zip
    """
}
