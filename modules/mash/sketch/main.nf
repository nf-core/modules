process MASH_SKETCH {
    tag "$meta.id"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::mash=2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1"
    } else {
        container "quay.io/biocontainers/mash:2.3--he348c14_1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.msh")        , emit: mash
    tuple val(meta), path("*.mash_stats") , emit: stats
    path "versions.yml"                   , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mash \\
        sketch \\
        $args \\
        -p $task.cpus \\
        -o ${prefix} \\
        -r $reads \\
        2> ${prefix}.mash_stats

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
