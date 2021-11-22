process GRAPHMAP2_INDEX {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0"
    } else {
        container "quay.io/biocontainers/graphmap:0.6.3--he513fc3_0"
    }

    input:
    path fasta

    output:
    path "*.gmidx"      , emit: index
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    graphmap2 \\
        align \\
        -t $task.cpus \\
        -I \\
        $args \\
        -r $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
