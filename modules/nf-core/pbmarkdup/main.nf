process PBMARKDUP {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pbmarkdup:1.2.0--h9ee0642_0' :
        'biocontainers/pbmarkdup:1.2.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: markduped
    tuple val(meta), path("${dupfile_name}")    , emit: dupfile   , optional: true
    tuple val(meta), path("*.pbmarkdup.log")    , emit: log       , optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args  ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    // To allow multiple input types/files: (compressed) fasta, fastq, bam; Determine suffix from input file names
    suffix        =
        input.find {
            it.name ==~ /.*\.(fasta|fa|fna)(\.gz)?$/ }?.with { f ->
            f.name.tokenize('.').takeRight(f.name.endsWith('.gz') ? 2 : 1).join('.')
        } ?:
        input.find { it.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }?.with { f ->
            f.name.tokenize('.').takeRight(f.name.endsWith('.gz') ? 2 : 1).join('.')
        } ?:
        input[0].extension
    dupfile_name = args.contains('--dup-file') ? (args =~ /--dup-file\s+(\S+)/)[0][1] : ''
    // PBmarkdup does not automatically gzip output files, even if the output file name ends with .gz.
    // Gzip the duplicate file in that case
    def compress_dup_args = (dupfile_name && dupfile_name.endsWith('.gz')) ?
        """
        if ! gzip -t "${dupfile_name}" 2>/dev/null; then
            gzip -c "${dupfile_name}" > tmp.dup.gz && mv "tmp.dup.gz" "${dupfile_name}"
        fi
        """ : ''
    def log_args = args.contains('--log-level') ? " > ${prefix}.pbmarkdup.log" : ''
    def file_list = input.collect { it.getName() }.join(' ')

    // Check file name collisions between input, output, and duplicate file
    if (file_list.contains("${prefix}.${suffix}"))
        error """Output file `${prefix}.${suffix}` conflicts with an input file.
        Please change the output `$prefix` or input file names."""
    if (dupfile_name) {
        if (file_list.contains(dupfile_name))
            error """Duplicate file `$dupfile_name` conflicts with an input file.
            Please change the duplicate file name `$dupfile_name` or input file names."""

        if (dupfile_name == "${prefix}.${suffix}")
            error """Duplicate file `$dupfile_name` cannot be the same as the output file name.
            Please change the duplicate file name `$dupfile_name` or output prefix `$prefix`."""
    }

    """
    pbmarkdup \\
        -j ${task.cpus} \\
        ${file_list} \\
        ${prefix}.${suffix} \\
        ${args} \\
        ${log_args}

    if [[ ${prefix}.${suffix} == *.gz ]] && ! gzip -t "${prefix}.${suffix}" 2>/dev/null; then
        gzip -c "${prefix}.${suffix}" > tmp.gz && mv "tmp.gz" "${prefix}.${suffix}"
    fi

    ${compress_dup_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmarkdup: \$(echo \$(pbmarkdup --version 2>&1) | awk 'BEFORE{FS=" "}{print \$2}')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args  ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    // To allow multiple input types/files: (compressed) fasta, fastq, bam; Determine suffix from input file names
    suffix        =
        input.find {
            it.name ==~ /.*\.(fasta|fa|fna)(\.gz)?$/ }?.with { f ->
            f.name.tokenize('.').takeRight(f.name.endsWith('.gz') ? 2 : 1).join('.')
        } ?:
        input.find { it.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }?.with { f ->
            f.name.tokenize('.').takeRight(f.name.endsWith('.gz') ? 2 : 1).join('.')
        } ?:
        input[0].extension
    dupfile_name  = args.contains('--dup-file') ? (args =~ /--dup-file\s+(\S+)/)[0][1] : ''
    def log_args  = args.contains('--log-level') ? " > ${prefix}.pbmarkdup.log" : ''
    def file_list = input.collect { it.getName() }.join(' ')
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmarkdup: \$(echo \$(pbmarkdup --version 2>&1) | awk 'BEFORE{FS=" "}{print \$2}')
    END_VERSIONS
    """
}
