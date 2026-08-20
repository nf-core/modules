process ODGI_LAYOUT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/odgi:0.9.0--py312h5e9d817_1':
        'quay.io/biocontainers/odgi:0.9.0--py312h5e9d817_1' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.lay"), optional: true, emit: lay
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    tuple val("${task.process}"), val('odgi'), eval("odgi version | sed 's/^v//; s/-.*//'"), emit: versions_odgi, topic: versions

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
    """

    stub:
    """
    touch ${meta.id}.lay
    """
}
