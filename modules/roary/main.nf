process ROARY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::roary=3.13.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0"
    } else {
        container "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0"
    }

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                    , emit: results
    tuple val(meta), path("results/*.aln"), optional: true, emit: aln
    path "versions.yml"                                   , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    roary \\
        $options.args \\
        -p $task.cpus \\
        -f results/ \\
        $gff

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( roary --version )
    END_VERSIONS
    """
}
