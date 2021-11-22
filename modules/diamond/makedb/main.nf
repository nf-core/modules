process DIAMOND_MAKEDB {
    tag "$fasta"
    label 'process_medium'

    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? 'bioconda::diamond=2.0.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container 'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0'
    }

    input:
    path fasta

    output:
    path "${fasta}.dmnd", emit: db
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    diamond \\
        makedb \\
        --threads $task.cpus \\
        --in  $fasta \\
        -d $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
