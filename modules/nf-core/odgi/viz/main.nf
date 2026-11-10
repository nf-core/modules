// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ODGI_VIZ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c519019f443966467de7ce3245a9bcc384e27b576cade1dafa675df199367c4/data'
:         'community.wave.seqera.io/library/odgi:0.9.0--48ccb29b570a8d65' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("*.png"), emit: png
    tuple val("${task.process}"), val('odgi'), eval("odgi version | sed 's/^v//; s/-.*//'"), emit: versions_odgi, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    odgi \\
        viz \\
        --threads $task.cpus \\
        --idx ${graph} \\
        --out ${prefix}.png \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png
    """
}
