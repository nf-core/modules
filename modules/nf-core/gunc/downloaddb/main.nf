process GUNC_DOWNLOADDB {
    tag "$db_name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/252ddb9a919c571e652a84faf3e07ea565f74e00307e5f1b24c26539f6020b61/data'
        : 'community.wave.seqera.io/library/gunc:1.1.0--98c05141813649a5'}"

    input:
    val db_name

    output:
    path "*.dmnd", emit: db
    tuple val("${task.process}"), val('gunc'), eval('gunc --version'), emit: versions_gunc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gunc download_db . -db $db_name $args
    """

    stub:

    """
    touch ${db_name}.dmnd
    """
}
