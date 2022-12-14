process MINIPROT_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::miniprot=0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniprot:0.5--h7132678_0':
        'quay.io/biocontainers/miniprot:0.5--h7132678_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mpi"), emit: index
    path "versions.yml"           , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniprot: \$(miniprot --version 2>&1)
    END_VERSIONS
    """
}
