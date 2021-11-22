process COOLER_ZOOMIFY {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0"
    }

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.mcool"), emit: mcool
    path "versions.yml"             , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cooler zoomify \\
        $options.args \\
        -n $task.cpus \\
        -o ${prefix}.mcool \\
        $cool

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
}
