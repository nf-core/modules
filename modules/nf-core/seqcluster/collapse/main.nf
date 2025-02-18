process SEQCLUSTER_COLLAPSE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqcluster:1.2.9--pyh5e36f6f_0':
        'biocontainers/seqcluster:1.2.9--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    path "versions.yml"                 , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcluster: \$(echo \$(seqcluster --version 2>&1) | sed 's/^.*seqcluster //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqcluster: \$(echo \$(seqcluster --version 2>&1) | sed 's/^.*seqcluster //')
    END_VERSIONS
    """
}
