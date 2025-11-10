process NANOFILT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0':
        'biocontainers/nanofilt:2.8.0--py_0' }"

    input:
    tuple val(meta), path(reads)
    path  summary_file

    output:
    tuple val(meta), path("*.fastq.gz"), emit: filtreads
    path "*.log"                       , optional: true, emit: log_file
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "filtered_${meta.id}"
    def sum    = summary_file ? "--summary ${summary_file}" : ''
    """
    gunzip \\
        -c $reads \\
        | NanoFilt \\
        $sum \\
        $args \\
        | gzip > ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$(NanoFilt -v | sed 's/NanoFilt //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "filtered_${meta.id}"
    def sum    = summary_file ? "--summary ${summary_file}" : ''
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$(NanoFilt -v | sed 's/NanoFilt //')
    END_VERSIONS
    """
}
