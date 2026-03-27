process GUNC_MERGECHECKM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/252ddb9a919c571e652a84faf3e07ea565f74e00307e5f1b24c26539f6020b61/data'
        : 'community.wave.seqera.io/library/gunc:1.1.0--98c05141813649a5'}"

    input:
    tuple val(meta), path(gunc_file), path(checkm_file)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('gunc'), eval('gunc --version'), emit: versions_gunc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gunc \\
        merge_checkm \\
        ${args} \\
        -g ${gunc_file} \\
        -c ${checkm_file} \\
        -o .
    """

    stub:
    """
    touch gunc_merge_checkm.tsv
    """
}
