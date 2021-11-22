process MINIMAP2_INDEX {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0"
    } else {
        container "quay.io/biocontainers/minimap2:2.21--h5bf99c6_0"
    }

    input:
    path fasta

    output:
    path "*.mmi"        , emit: index
    path "versions.yml" , emit: versions

    script:
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
