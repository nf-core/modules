process PICARD_FASTQTOSAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
        'biocontainers/picard:3.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        log.warn '[Picard FastqToSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    }
    def avail_mem = task.memory ? (task.memory.mega*0.8).intValue() : 3072
    def input = meta.single_end ? "--FASTQ ${reads}" : "--FASTQ ${reads[0]} --FASTQ2 ${reads[1]}"
    def sample_name = args.contains("--SAMPLE_NAME") || args.contains("-SM") ? "" : "--SAMPLE_NAME ${prefix}"
    """
    picard \\
        -Xmx${avail_mem}M \\
        FastqToSam \\
        ${args} \\
        ${input} \\
        ${sample_name} \\
        --OUTPUT ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard FastqToSam --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
