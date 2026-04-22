process AGAT_SPFILTERBYORFSIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("*.passed.gff"), emit: passed_gff
    tuple val(meta), path("*.failed.gff"), emit: failed_gff
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def config_arg  = config ? "-c ${config}" : ''
    if( "${gxf}" in [ "${prefix}.passed.gff", "${prefix}.failed.gff" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    agat_sp_filter_by_ORF_size.pl \\
        -g ${gxf} \\
        ${args} \\
        ${config_arg} \\
        -o ${prefix}

    mv \\
        ${prefix}_NOT* \\
        "${prefix}.failed.gff"

    mv \\
        ${prefix}_* \\
        "${prefix}.passed.gff"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "${gxf}" in [ "${prefix}.passed.gff", "${prefix}.failed.gff" ] ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.passed.gff
    touch ${prefix}.failed.gff
    """
}
