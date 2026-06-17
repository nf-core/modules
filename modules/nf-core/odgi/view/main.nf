process ODGI_VIEW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.9.0--py312h5e9d817_1':
        'quay.io/biocontainers/odgi:0.9.0--py312h5e9d817_1' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.gfa"), emit: gfa
    tuple val("${task.process}"), val('odgi'), eval("odgi version | sed 's/^v//; s/-.*//'"), emit: versions_odgi, topic: versions


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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gfa
    """
}
