process FOLDDISCO_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/folddisco:2.9375a2d--hb42e459_0':
        'quay.io/biocontainers/folddisco:2.9375a2d--hb42e459_0' }"

    input:
    tuple val(meta), path(structures, stageAs: "structures/*")

    output:
    tuple val(meta), path("${prefix}"), emit: index
    tuple val("${task.process}"), val('folddisco'), eval("folddisco version"), topic: versions, emit: versions_folddisco

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    folddisco \\
        index \\
        --pdbs structures \\
        --index ${prefix}/${prefix} \\
        --threads $task.cpus \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}
    touch ${prefix}/${prefix}.offset
    touch ${prefix}/${prefix}.lookup
    touch ${prefix}/${prefix}.type
    """
}
