process TSEBRA {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tsebra:1.1.2.5--pyhca03a8a_0':
        'quay.io/biocontainers/tsebra:1.1.2.5--pyhca03a8a_0' }"

    input:
    tuple val(meta), path(gtfs)
    path hints_files
    path keep_gtfs
    path config

    output:
    tuple val(meta), path("*.gtf"), emit: tsebra_gtf
    tuple val(meta), path("*.tsv"), emit: tsebra_scores
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('tsebra'), val("1.1.2.5"), emit: versions_tsebra, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args                                     ?: ''
    def prefix      = task.ext.prefix                                   ?: "${meta.id}"
    def gtf_arg     = '-g ' + gtfs.collect { gtf -> "$gtf" }.join(',')
    def hints_arg   = '-e ' + hints_files.collect { hint -> "$hint" }.join(',')
    def keep_arg    = keep_gtfs                                         ? ( '-k ' + keep_gtfs.collect { gtf -> "$gtf" }.join(',') ) : ''
    def config_arg  = config                                            ? "-c $config"                                      : ''
    """
    tsebra.py \\
        ${gtf_arg} \\
        ${hints_arg} \\
        ${keep_arg} \\
        ${config_arg} \\
        ${args} \\
        -o ${prefix}.gtf \\
        -s ${prefix}.tsv
    """

    stub:
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    touch ${prefix}.tsv
    """
}
