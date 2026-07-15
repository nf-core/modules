process TREERECS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/treerecs:1.2--h9f5acd7_3':
        'quay.io/biocontainers/treerecs:1.2--h9f5acd7_3' }"

    input:
    tuple val(meta), path(genetree), path(speciestree)

    output:
    tuple val(meta), path("*.xml"), emit: rectree
    tuple val("${task.process}"), val('treerecs'), eval("treerecs --version | sed -n 's/^Treerecs \\([0-9.]*\\)\$/\\1/p'"), topic: versions, emit: versions_treerecs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    treerecs \\
        $args \\
        --genetree $genetree \\
        --speciestree $speciestree \\
        --outdir treerecs_tmp \\
        --output-format recphyloxml \\
        --force \\
        --parallelize

    mv treerecs_tmp/* ${prefix}.recphylo.xml
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.recphylo.xml
    """
}
