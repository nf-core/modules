process PIRATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pirate=1.0.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate%3A1.0.4--hdfd78af_1' :
        'quay.io/biocontainers/pirate:1.0.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                   , emit: results
    tuple val(meta), path("results/core_alignment.fasta"), optional: true, emit: aln
    path "versions.yml"                                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    PIRATE \\
        $args \\
        --threads $task.cpus \\
        --input ./ \\
        --output results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
