// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CENTRIFUGE_KREPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/273bd9423f536f18a99e522fd2870d257f0f39d6b24329a8e0356c95e3ab1e95/data'
        : 'community.wave.seqera.io/library/centrifuge:1.0.4.2--d96245be9f9eb1f0'}"

    input:
    tuple val(meta), path(report)
    path db

    output:
    tuple val(meta), path('*.txt'), emit: kreport
    tuple val("${task.process}"), val("centrifuge"), eval("centrifuge --version 2>&1 | sed '1!d;s/.* version //'"), emit: versions_centrifuge, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_name=`find -L ${db} -name "*.1.cf" -not -name "._*"  | sed 's/\\.1.cf\$//'`
    centrifuge-kreport \\
        ${args} \\
        -x \$db_name \\
        ${report} > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
