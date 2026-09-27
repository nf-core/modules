process PBTK_BAM2FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbtk:3.1.1--h9ee0642_0':
        'quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(pbi)

    output:
    tuple val(meta), path("*.$extension")   , emit: fastq
    tuple val("${task.process}"), val('pbtk'), eval("bam2fastq --version 2>&1 | sed -n 's/.*bam2fastq //p' | cut -d' ' -f1"), emit: versions_bam2fastq, topic: versions

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
    """

    stub:
    def args    = task.ext.args         ?: ''
    def prefix  = task.ext.prefix       ?: "${meta.id}"
    def gzip    = args.contains('-u')   ? ''            : "gzip *.fastq"
    extension   = args.contains('-u')   ? 'fastq'       : 'fastq.gz'
    """
    touch ${prefix}.fastq
    $gzip
    """
}
