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
    suffix       = input[0].getExtension()             // To allow multiple input types
    dupfile_name = args.contains('--dup-file') ? (args =~ /--dup-file\s+(\S+)/)[0][1] : ''
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
        $args \\
        ${log_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmarkdup: \$(echo \$(pbmarkdup --version 2>&1) | awk 'BEFORE{FS=" "}{print \$2}')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args  ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    suffix        = input[0].getExtension()             // To allow multiple input types
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
