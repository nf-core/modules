// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process INSITUGAME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insitugame:1.0.7--pyhdfd78af_0':
        'quay.io/biocontainers/insitugame:1.0.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(spatial_data)

    output:
    tuple val(meta), path("*.png") , emit: figures, optional: true
    tuple val(meta), path("*.json"), emit: results
    tuple val(meta), path("*.csv") , emit: tables
    tuple val("${task.process}"), val('insitugame'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('insitugame'))\""), topic: versions, emit: versions_insitugame

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    insitugame \\
        $args \\
        --input $spatial_data \\
        --output-prefix ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.json
    touch ${prefix}_tables.csv
    touch ${prefix}_figure.png
    """
}
