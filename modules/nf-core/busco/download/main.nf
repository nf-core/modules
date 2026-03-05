process BUSCO_DOWNLOAD {
    tag "$lineage"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c607f319867d96a38c8502f751458aa78bbd18fe4c7c4fa6b9d8350e6ba11ebe/data'
        : 'community.wave.seqera.io/library/busco_sepp:f2dbc18a2f7a5b64'}"

    input:
    val lineage

    output:
    path "busco_downloads", emit: download_dir
    tuple val("${task.process}"), val('busco'), eval("busco --version 2> /dev/null | sed 's/BUSCO //g'"), emit: versions_busco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    busco \\
        --download $lineage \\
        $args
    """

    stub:
    """
    mkdir busco_downloads
    """
}
