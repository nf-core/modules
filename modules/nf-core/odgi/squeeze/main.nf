// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ODGI_SQUEEZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c519019f443966467de7ce3245a9bcc384e27b576cade1dafa675df199367c4/data'
:         'community.wave.seqera.io/library/odgi:0.9.0--48ccb29b570a8d65' }"

    input:
    tuple val(meta), path(graphs)

    output:
    tuple val(meta), path("*.og"), emit: graph
    tuple val("${task.process}"), val('odgi'), eval("odgi version | sed 's/^v//; s/-.*//'"), emit: versions_odgi, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.og
    """
}
