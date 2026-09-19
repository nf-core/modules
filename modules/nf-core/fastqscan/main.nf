process FASTQSCAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-scan:0.4.4--h7d875b9_0' :
        'quay.io/biocontainers/fastq-scan:0.4.4--h7d875b9_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val('fastqscan'), eval('fastq-scan -v 2>&1 | sed \'s/^.*fastq-scan //\''), emit: versions_fastqscan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $reads | \\
        fastq-scan \\
        $args > ${prefix}.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json
    """
}
