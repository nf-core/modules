process GRAPHMAP2_INDEX {
    label 'process_medium'

    conda "bioconda::graphmap=0.6.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphmap:0.6.3--he513fc3_0' :
        'biocontainers/graphmap:0.6.3--he513fc3_0' }"

    input:
    path fasta

    output:
    path "*.gmidx"      , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
    "${task.process}":
        graphmap2: \$(echo \$(graphmap2 align 2>&1) | sed 's/^.*Version: v//; s/ .*\$//')
    END_VERSIONS
    """
}
