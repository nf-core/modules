process PICARD_SPLITSAMBYNUMBEROFREADS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/picard%3A3.4.0--hdfd78af_0'
        : 'biocontainers/picard:3.4.0--hdfd78af_0'}"

    input:
    tuple val(meta), path(bam)
    val split_to_N_reads
    val split_to_N_files
    file arguments_file

    output:
    tuple val(meta), path("picardsplit/*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_file = arguments_file ? "--arguments_file ${arguments_file}" : ""
    assert (!split_to_N_reads && split_to_N_files) || (split_to_N_reads && !split_to_N_files) : "You must provide either 'split_to_N_reads' or 'split_to_N_files', but not both."
    def split_arg = ''
    if (!split_to_N_reads) {
        split_arg = "--SPLIT_TO_N_FILES ${split_to_N_files}"
    }
    else {
        split_arg = "--SPLIT_TO_N_READS ${split_to_N_reads}"
    }
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard SplitSamByNumberOfReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """

    mkdir picardsplit

    picard \\
        -Xmx${avail_mem}M \\
        SplitSamByNumberOfReads \\
        ${split_arg} \\
        ${args_file} \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT picardsplit \\
        --OUT_PREFIX ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SplitSamByNumberOfReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:n)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_file = arguments_file ? "--arguments_file ${arguments_file}" : ""
    assert (!split_to_N_reads && split_to_N_files) || (split_to_N_reads && !split_to_N_files) : "You must provide either 'split_to_N_reads' or 'split_to_N_files', but not both."
    def split_arg = ''
    if (!split_to_N_reads) {
        split_arg = "--SPLIT_TO_N_FILES ${split_to_N_files}"
    }
    else {
        split_arg = "--SPLIT_TO_N_READS ${split_to_N_reads}"
    }
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[Picard SplitSamByNumberOfReads] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    mkdir picardsplit
    touch picardsplit/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard SplitSamByNumberOfReads --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
