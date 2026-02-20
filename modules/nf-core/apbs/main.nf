process APBS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/288b945c2895993a166b5baf67efd85f2111a3efe64bd189dae731e697eae9ef/data':
        'community.wave.seqera.io/library/apbs:3.4.1--298b75172827aae7' }"

    input:
    tuple val(meta), path(input_file) // TODO add the pqr file, path (in), path (pqr)

    output:
    tuple val(meta), path("*.dx") , emit: dx
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('apbs'), eval("apbs --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | tail -n1"), topic: versions, emit: versions_apbs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    apbs \\
        $args \\
        ${input_file} \\
        2>&1 | tee ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args} >/dev/null

    touch ${prefix}.dx
    touch ${prefix}.log
    """
}
