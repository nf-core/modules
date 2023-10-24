process ODGI_VIEW {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::odgi=0.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.8.3--py310h6cc9453_0':
        'biocontainers/odgi:0.8.3--py310h6cc9453_0' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi \\
        view \\
        --threads $task.cpus \\
        --idx ${graph} \\
        --to-gfa \\
        $args > ${prefix}.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(echo \$(odgi version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
