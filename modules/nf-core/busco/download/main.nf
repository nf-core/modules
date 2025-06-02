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
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${lineage}"
    """
    busco \\
        --download $lineage \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2> /dev/null | sed 's/BUSCO //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir busco_downloads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2> /dev/null | sed 's/BUSCO //g' )
    END_VERSIONS
    """
}
