process MINIPROT_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniprot:0.11--he4a0461_2':
        'biocontainers/miniprot:0.11--he4a0461_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mpi"), emit: index
    tuple val("${task.process}"), val('miniprot'), eval("miniprot --version"), topic: versions, emit: versions_miniprot

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    miniprot \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mpi \\
        $args \\
        $fasta
    """

    stub:
    """
    touch ${fasta.baseName}.mpi
    """
}
