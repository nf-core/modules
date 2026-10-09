process CHOPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b21cace9a3bd19a10af876ac9f628117d1a16264ded6b76c3da0b762d2e78cb5/data':
        'community.wave.seqera.io/library/chopper_pigz:26723dce47edae8f' }"

    input:
    tuple val(meta), path(fastq)
    path  fasta

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    tuple val("${task.process}"), val('chopper'), eval("chopper --version 2>&1 | cut -d ' ' -f 2"), topic: versions, emit: versions_chopper
    tuple val("${task.process}"), val('pigz'),    eval('pigz --version | sed "s/pigz //g"'),        topic: versions, emit: versions_chopper

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def args3  = task.ext.args3  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_filtering = fasta ? "--contam ${fasta}" : ""

    if ("$fastq" == "${prefix}.fastq.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    pigz \\
        $args \\
        $fastq | \\
    chopper \\
        --threads $task.cpus \\
        $fasta_filtering \\
        $args2 | \\
    pigz \\
        $args3 > ${prefix}.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | pigz > ${prefix}.fastq.gz
    """
}
