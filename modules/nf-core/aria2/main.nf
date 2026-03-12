process ARIA2 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c0d3d867f5bc805b926b08ee761a993b24062739743eb82cc56363e0f7817d/data' :
        'community.wave.seqera.io/library/aria2:1.37.0--3a9ec328469995dd' }"

    input:
    tuple val(meta), val(source_url)

    output:
    tuple val(meta), path("$downloaded_file"), emit: downloaded_file
    tuple val("${task.process}"), val("aria2"), eval("aria2c --version 2>&1 | sed -n 's/^aria2 version \\([^ ]*\\).*/\\1/p'"), emit: versions_aria2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    downloaded_file = source_url.split("/")[-1]

    """
    aria2c \\
        --check-certificate=false \\
        ${args} \\
        ${source_url}
    """

    stub:
    def args = task.ext.args ?: ''
    downloaded_file = source_url.split("/")[-1]

    """
    echo ${args}

    touch ${downloaded_file}
    """
}
