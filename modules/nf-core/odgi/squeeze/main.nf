process ODGI_SQUEEZE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::odgi=0.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.8.3--py310h6cc9453_0':
        'biocontainers/odgi:0.8.3--py310h6cc9453_0' }"

    input:
    tuple val(meta), path(graphs)

    output:
    tuple val(meta), path("*.og"), emit: graph
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls *.og > files
    odgi \\
        squeeze \\
        $args \\
        --threads $task.cpus \\
        --input-graphs files \\
        -o ${prefix}.og


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(echo \$(odgi version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v'))
    END_VERSIONS
    """
}
