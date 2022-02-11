process MINIMAP2_INDEX {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0' :
        'quay.io/biocontainers/minimap2:2.21--h5bf99c6_0' }"

    input:
    path fasta

    output:
    path "*.mmi"        , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
