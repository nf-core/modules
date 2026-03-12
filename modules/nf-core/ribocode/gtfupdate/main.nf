process RIBOCODE_GTFUPDATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe815db0864b45b91afc7bc84c55cb60acb0035e7248dda7f480a55c4cb105d7/data':
        'community.wave.seqera.io/library/ribocode:1.2.15--5530b252f5433a62' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.gtf")                                                  , emit: gtf
    tuple val("${task.process}"), val('ribocode'), eval('RiboCode --version  2>&1') , emit: versions_ribocode, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    GTFupdate \\
        ${gtf} \\
        $args \\
        > ${prefix}_updated.gtf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_updated.gtf
    """
}
