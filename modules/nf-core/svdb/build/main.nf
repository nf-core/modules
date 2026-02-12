process SVDB_BUILD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ff/ff995c756aa8a3c0af13b1d054eacd536a11d35de5fa288dacf558bc21696968/data':
        'community.wave.seqera.io/library/bcftools_svdb:ae3b14d2d608fd81' }"

    input:
    tuple val(meta), path(input)
    val(input_type)

    output:
    tuple val(meta), path("*.db"), emit: db
    tuple val("${task.process}"), val('svdb'), eval("svdb | sed -nE 's/.*SVDB-([0-9.]+).*/\\1/p'"), emit: versions_svdb, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    if (!input_type.matches('folder|files')) { error "Unrecognised input type. Options are: 'folder', 'files'" }

    """
    svdb \\
        --build \\
        ${args} \\
        --$input_type $input \\
        --prefix ${prefix}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.db

    """
}
