process SEQCLUSTER_COLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqcluster:1.2.9--pyh5e36f6f_0':
        'quay.io/biocontainers/seqcluster:1.2.9--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    tuple val("${task.process}"), val('seqcluster'), eval("seqcluster --version 2>&1 | tail -n 1 | sed 's/^seqcluster //'"), topic: versions, emit: versions_seqcluster

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fastq" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    seqcluster \\
        collapse \\
        $args \\
        -f $fastq  \\
        -o collapsed

    gzip collapsed/*_trimmed.fastq
    mv collapsed/*_trimmed.fastq.gz ${prefix}.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
