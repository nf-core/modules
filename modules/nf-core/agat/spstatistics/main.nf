process AGAT_SPSTATISTICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.6.1--pl5321hdfd78af_1' :
        'biocontainers/agat:1.6.1--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.txt"), emit: stats_txt
    tuple val("${task.process}"), val('agat'), eval("agat --version | sed 's/v//'"), topic: versions, emit: versions_agat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_sp_statistics.pl \\
        --gff ${gff} \\
        --output ${prefix}.stats.txt \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats.txt
    """
}
