process PBTK_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0':
        'biocontainers/pbtk:3.1.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(pbi)

    output:
    tuple val(meta), path("*.$extension")   , emit: fastq
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args         ?: ''
    def prefix  = task.ext.prefix       ?: "${meta.id}"
    extension   = args.contains('-u')   ? 'fastq'       : 'fastq.gz'
    """
    bam2fastq \\
        $args \\
        -j $task.cpus \\
        -o ${prefix} \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq: \$(bam2fastq --version 2>&1 | sed -n 's/.*bam2fastq \\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args         ?: ''
    def prefix  = task.ext.prefix       ?: "${meta.id}"
    def gzip    = args.contains('-u')   ? ''            : "gzip *.fastq"
    extension   = args.contains('-u')   ? 'fastq'       : 'fastq.gz'
    """
    touch ${prefix}.fastq
    $gzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam2fastq: \$(bam2fastq --version 2>&1 | sed -n 's/.*bam2fastq \\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
