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
    path "versions.yml"                   , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
