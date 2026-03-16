process AGAT_SPKEEPLONGESTISOFORM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gxf)
    path config

    output:
    tuple val(meta), path("${output}"), emit: gff
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ""
    output           = "${prefix}.longest.gff"
    """
    agat_sp_keep_longest_isoform.pl \\
        --gff ${gxf} \\
        ${config_param} \\
        --out ${output} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output     = "${prefix}.longest.gff"
    """
    touch ${output}
    """
}
