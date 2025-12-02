
process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.8--py39h2de1943_0':
        'biocontainers/whatshap:2.8--py39h2de1943_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv"),                                            emit: tsv
    tuple val(meta), path("*.txt"),                                            emit: txt
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"), emit: versions_whatshap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    whatshap stats \\
        $args \\
        --sample ${meta.id} \\
        --tsv ${prefix}.tsv \\
        $vcf \\
        > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    touch ${prefix}.txt
    """
}
