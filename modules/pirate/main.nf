process PIRATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pirate=1.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pirate%3A1.0.4--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/pirate:1.0.4--hdfd78af_1"
    }

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                   , emit: results
    tuple val(meta), path("results/core_alignment.fasta"), optional: true, emit: aln
    path "versions.yml"                                                  , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    PIRATE \\
        $options.args \\
        --threads $task.cpus \\
        --input ./ \\
        --output results/

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
