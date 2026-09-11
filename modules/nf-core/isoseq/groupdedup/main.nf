process ISOSEQ_GROUPDEDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:26.2.0--h9ee0642_0' :
        'quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.dedup.bam")     , emit: bam
    tuple val(meta), path("*.dedup.bam.pbi") , emit: pbi
    tuple val(meta), path("*.dedup.fasta.gz"), emit: fasta
    tuple val("${task.process}"), val('isoseq'), eval("isoseq groupdedup --version | sed -n '1s/isoseq groupdedup \\([0-9.]*\\).*/\\1/p'"), emit: versions_isoseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isoseq \\
        groupdedup \\
        -j $task.cpus \\
        $args \\
        $bam \\
        ${prefix}.dedup.bam

    gzip ${prefix}.dedup.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dedup.bam
    touch ${prefix}.dedup.bam.pbi
    echo "" | gzip > ${prefix}.dedup.fasta.gz
    """
}
