process SVDSS_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdss:2.1.1--he17396a_0' :
        'quay.io/biocontainers/svdss:2.1.1--he17396a_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.fmd"), emit: fmd
    tuple val("${task.process}"), val('svdss'), eval("SVDSS --version 2>&1 | sed 's/SVDSS, //'"), emit: versions_svdss, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    SVDSS index \\
        -t ${task.cpus} \\
        -d ${fasta} \\
        ${args} \\
        > ${prefix}.fmd
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fmd
    """
}
