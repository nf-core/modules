process GRAPHMAP2_INDEX {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::graphmap=0.6.3" : null)
    def container_image = "/graphmap:0.6.3--he513fc3_0"
                                                  container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

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
