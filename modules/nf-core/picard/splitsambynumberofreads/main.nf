process PICARD_SPLITSAMBYNUMBEROFREADS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0' :
        'biocontainers/picard:3.4.0--hdfd78af_0'}"

    input:
    tuple val(meta) , path(input)
    tuple val(meta2), path(fasta), path(fai)
    val split_to_N_reads
    val split_to_N_files
    path arguments_file

    output:
    tuple val(meta), path("picardsplit/*.{bam,sam,cram}"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
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
        ${reference} \\
        ${args_file} \\
        ${args} \\
        --INPUT ${input} \\
        --OUTPUT picardsplit \\
        --OUT_PREFIX ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard PositionBasedDownsampleSam --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!split_to_N_reads) {
        split_arg = "--SPLIT_TO_N_FILES ${split_to_N_files}"
    }
    else {
        split_arg = "--SPLIT_TO_N_READS ${split_to_N_reads}"
    }
    def avail_mem = 3072
    def ext = input.name.tokenize('.')[-1]
    if ( !['bam', 'cram', 'sam'].contains(ext) )
        ext = 'bam'
    """
    mkdir picardsplit

    if echo "${split_to_N_files}" | grep -qE '^[1-9][0-9]*\$'; then
        i=1
        while [ \$i -le ${split_to_N_files} ]; do
            fname=\$(printf "%04d" \$i)
            touch picardsplit/${prefix}_\${fname}.${ext}
            i=\$((i + 1))
        done
    else
        touch picardsplit/${prefix}_0001.${ext}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard PositionBasedDownsampleSam --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
