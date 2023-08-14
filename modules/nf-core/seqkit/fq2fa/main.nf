process SEQKIT_FQ2FA {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqkit=2.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.5.0--h9ee0642_0':
        'biocontainers/seqkit:2.5.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fa.gz"), emit: fasta
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        fq2fa \\
        $args \\
        -j $task.cpus \\
        -o ${prefix}.fa.gz \\
        $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(seqkit version 2>&1) | sed 's/seqkit v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(seqkit version 2>&1) | sed 's/seqkit v//' ))
    END_VERSIONS
    """
}
