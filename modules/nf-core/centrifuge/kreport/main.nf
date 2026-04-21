process CENTRIFUGE_KREPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4.2--hdcf5f25_0'
        : 'quay.io/biocontainers/centrifuge:1.0.4.2--hdcf5f25_0'}"

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
