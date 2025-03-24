process FASTQDL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-dl:3.0.0--pyhdfd78af_0':
        'biocontainers/fastq-dl:3.0.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("test"), emit: db
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-dl \\
        $args \\
        --accession $accession \\
        --cpus $task.cpus \\
        --outdir ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version |& sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir ${prefix}
    echo "" | gzip > ${accession}.fastq.gz
    echo "" | gzip > ${accession}_1.fastq.gz
    echo "" | gzip > ${accession}_2.fastq.gz
    touch ${prefix}-run-info.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-dl: \$(fastq-dl --version |& sed 's/.* //')
    END_VERSIONS
    """
}
