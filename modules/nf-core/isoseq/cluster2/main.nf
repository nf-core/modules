process ISOSEQ_CLUSTER2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoseq:26.2.0--h9ee0642_0' :
        'quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.transcripts.bam")               , emit: bam
    tuple val(meta), path("*.transcripts.bam.pbi")           , emit: pbi
    tuple val(meta), path("*.transcripts.cluster_report.csv"), emit: cluster_report
    tuple val("${task.process}"), val('isoseq'), eval("isoseq cluster2 --version | sed -n '1s/isoseq cluster2 \\([0-9.]*\\).*/\\1/p'"), emit: versions_isoseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    isoseq \\
        cluster2 \\
        -j $task.cpus \\
        $args \\
        $bam \\
        ${prefix}.transcripts.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.bam
    touch ${prefix}.transcripts.bam.pbi
    touch ${prefix}.transcripts.cluster_report.csv
    """
}
