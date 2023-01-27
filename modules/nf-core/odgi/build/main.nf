process ODGI_BUILD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::odgi=0.8.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.8.2--py310hc8f18ef_0':
        'quay.io/biocontainers/odgi:0.8.2--py310hc8f18ef_0' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.og"), emit: og
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi \\
        build \\
        --threads $task.cpus \\
        --gfa ${graph} \\
        --out ${prefix}.og \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(echo \$(odgi version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
