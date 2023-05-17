process ODGI_LAYOUT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::odgi=0.8.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.8.2--py310hc8f18ef_0':
        'quay.io/biocontainers/odgi:0.8.2--py310hc8f18ef_0' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.lay"), optional: true, emit: lay
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "--out ${meta.id}.lay"
    """
    odgi \\
        layout \\
        --threads $task.cpus \\
        --idx ${graph} \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(echo \$(odgi version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
