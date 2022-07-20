process PICARD_FASTQTOSAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picard=2.27.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.27.4--hdfd78af_0' }"

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
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard FastqToSam] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def input = meta.single_end ? "--FASTQ ${reads}" : "--FASTQ ${reads[0]} --FASTQ2 ${reads[1]}"
    def sample_name = args.contains("--SAMPLE_NAME") || args.contains("-SM") ? "" : "--SAMPLE_NAME ${prefix}"
    """
    picard \\
        -Xmx${avail_mem}g \\
        FastqToSam  \\
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
