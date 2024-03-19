process NANOFILT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0':
        'biocontainers/nanofilt:2.8.0--py_0' }"

    input:
    tuple val(meta), path(reads)
    val   readlength
    val   readqual

    output:
    tuple val(meta), path("*.fastq.gz"), emit: filtreads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "filtered_${meta.id}"
    """
    gunzip \\
        -c $reads \\
        | NanoFilt \\
        --maxlength $readlength \\
        -q $readqual \\
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
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanofilt: \$(NanoFilt -v | sed 's/NanoFilt //')
    END_VERSIONS
    """
}
