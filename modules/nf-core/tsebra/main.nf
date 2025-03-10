process TSEBRA {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tsebra:1.1.2.5--pyhca03a8a_0':
        'biocontainers/tsebra:1.1.2.5--pyhca03a8a_0' }"

    input:
    tuple val(meta), path(gtfs)
    path hints_files
    path keep_gtfs
    path config

    output:
    tuple val(meta), path("*.gtf"), emit: tsebra_gtf
    tuple val(meta), path("*.tsv"), emit: tsebra_scores
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args                                     ?: ''
    def prefix      = task.ext.prefix                                   ?: "${meta.id}"
    def gtf_arg     = '-g ' + gtfs.collect { "$it" }.join(',')
    def hints_arg   = '-e ' + hints_files.collect { "$it" }.join(',')
    def keep_arg    = keep_gtfs                                         ? ( '-k ' + keep_gtfs.collect { "$it" }.join(',') ) : ''
    def config_arg  = config                                            ? "-c $config"                                      : ''
    def VERSION     = '1.1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tsebra.py \\
        $gtf_arg \\
        $hints_arg \\
        $keep_arg \\
        $config_arg \\
        $args \\
        -o ${prefix}.gtf \\
        -s ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tsebra: $VERSION
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    def VERSION     = '1.1.2.5' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.gtf
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tsebra: $VERSION
    END_VERSIONS
    """
}
