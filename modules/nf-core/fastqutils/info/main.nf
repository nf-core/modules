process FASTQUTILS_INFO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq_utils:0.25.2--h96c455f_2':
        'biocontainers/fastq_utils:0.25.2--h96c455f_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.txt"), emit: info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq_info \\
        $args \\
        ${reads}

    echo "fastq_utils fastq_info ran and found no issues with ${reads}" > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_utils: \$(fastq_info -h 2>&1 | head -n 1 | sed 's/^fastq_utils //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq_utils: \$(fastq_info -h 2>&1 | head -n 1 | sed 's/^fastq_utils //')
    END_VERSIONS
    """
}
