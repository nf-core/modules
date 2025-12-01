
process WHATSHAP_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatshap:2.8--py39h2de1943_0':
        'quay.io/biocontainers/whatshap:2.8--py39h2de1943_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.stats.tsv"),                                                                         emit: stats
    tuple val(meta), path("*.blocks.gtf.gz"),                                                                     emit: blocks
    tuple val(meta), path("*.blocks.gtf.gz.tbi"),                                                                 emit: blocks_index
    tuple val("${task.process}"), val('whatshap'), eval("whatshap --version"),                                    emit: versions_whatshap, topic: versions
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | head -n 1 | sed 's/bgzip (htslib) //g'"), emit: versions_bgzip,    topic: versions
    tuple val("${task.process}"), val('tabix'), eval("tabix --version | head -n 1 | sed 's/tabix (htslib) //g'"), emit: versions_tabix,    topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
 
    """
    whatshap stats \\
        $args \\
        --sample ${meta.id} \\
        --tsv ${prefix}.stats.tsv \\
        --gtf ${prefix}.blocks.gtf \\
        $vcf

    bgzip \\
        -@ $task.cpus \\
        ${prefix}.blocks.gtf

    tabix \\
        ${prefix}.blocks.gtf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.stats.tsv
    echo | gzip > ${prefix}.blocks.gtf.gz
    touch ${prefix}.blocks.gtf.gz.tbi
    """
}
