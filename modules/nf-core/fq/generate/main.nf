process FQ_GENERATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.9.1--h9ee0642_0':
        'biocontainers/fq:0.9.1--h9ee0642_0' }"

    input:
    val meta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fq generate \\
        $args \\
        ${prefix}_R1.fastq.gz ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq generate --version | sed 's/fq-generate //g'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}_R1.fastq.gz
    echo | gzip > ${prefix}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq generate --version | sed 's/fq-generate //g'))
    END_VERSIONS
    """
}
