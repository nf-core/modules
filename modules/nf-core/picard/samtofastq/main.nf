process PICARD_SAMTOFASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.1.0--hdfd78af_0' :
        'biocontainers/picard:3.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq"), emit: reads
    tuple val("${task.process}"), val('picard'), eval("picard SamToFastq --version 2>&1 | sed -n 's/.*Version://p'"), topic: versions, emit: versions_picard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = meta.single_end ? "--FASTQ ${prefix}.fastq.gz" : "--FASTQ ${prefix}_1.fastq.gz --SECOND_END_FASTQ ${prefix}_2.fastq.gz"
    if (!task.memory) {
        log.warn '[Picard SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    """
    picard \\
        -Xmx${avail_mem}M \\
        SamToFastq \\
        ${args} \\
        --INPUT ${bam} \\
        ${output}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = meta.single_end ? "--FASTQ ${prefix}.fastq.gz" : "--FASTQ ${prefix}_1.fastq.gz --SECOND_END_FASTQ ${prefix}_2.fastq.gz"
    if (!task.memory) {
        log.warn '[Picard SamToFastq] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    """
    echo "picard \\
        -Xmx${avail_mem}M \\
        SamToFastq \\
        ${args} \\
        --INPUT ${bam} \\
        ${output}"
    touch ${output}
    """
}
