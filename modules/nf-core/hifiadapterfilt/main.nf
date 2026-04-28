process HIFIADAPTERFILT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiadapterfilt:3.0.0--hdfd78af_0':
        'quay.io/biocontainers/hifiadapterfilt:3.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.filt.fastq.gz")       , emit: fastq
    tuple val(meta), path("${prefix}.stats")               , emit: stats
    tuple val(meta), path("${prefix}.contaminant.blastout"), emit: blastout
    tuple val(meta), path("${prefix}.blocklist")           , emit: blocklist
    tuple val("${task.process}"), val('hifiadapterfilt'), eval("hifiadapterfilt.sh --version 2>&1 | head -1"), topic: versions, emit: versions_hifiadapterfilt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    ln -s \$(realpath ${bam}) ${prefix}.bam

    hifiadapterfilt.sh \\
        -t ${task.cpus} \\
        -o . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.filt.fastq.gz
    touch ${prefix}.stats
    touch ${prefix}.contaminant.blastout
    touch ${prefix}.blocklist
    """
}
