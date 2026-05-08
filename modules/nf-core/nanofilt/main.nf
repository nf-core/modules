process NANOFILT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0':
        'quay.io/biocontainers/nanofilt:2.8.0--py_0' }"

    input:
    tuple val(meta), path(reads)
    path  summary_file

    output:
    tuple val(meta), path("*.fastq.gz"), emit: filtreads
    path "*.log"                       , optional: true, emit: log_file
    tuple val("${task.process}"), val('nanofilt'), eval('NanoFilt -v | sed "s/NanoFilt //"'), emit: versions_nanofilt, topic: versions

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

    """

    stub:
    def prefix = task.ext.prefix ?: "filtered_${meta.id}"
    """
    echo "" | gzip > ${prefix}.fastq.gz
    touch ${prefix}.log

    """
}
