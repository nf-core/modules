process AGRVATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agrvate:1.0.2--hdfd78af_0' :
        'biocontainers/agrvate:1.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}-results/${fasta.baseName}-summary.tab"), emit: summary
    path "${fasta.baseName}-results"                                                , emit: results_dir
    tuple val("${task.process}"), val('agrvate'), eval("agrvate --version | sed 's/[^0-9.]//g'"), emit: versions_agrvate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    """
    agrvate \\
        ${args} \\
        -i ${fasta}
    """

    stub:

    """
    mkdir ${fasta.baseName}-results
    touch ${fasta.baseName}-results/${fasta.baseName}-summary.tab
    """
}
