// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CRAMINO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a1e9861d7c3ac6fd4290bb369746f68fa3cf773196be01e986f59228892d8d16/data'
        : 'community.wave.seqera.io/library/cramino:2.0.0--9dc31ae276d5fec3'}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt")  , emit: stats
    tuple val(meta), path("*.arrow"), optional: true, emit: arrow
    tuple val("${task.process}"), val('cramino'), eval("cramino -V | sed 's/cramino //'"), emit: versions_cramino, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cramino \\
        $args \\
        --threads $task.cpus \\
        ${bam} > ${prefix}.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arrow_file = args.contains('--arrow') ? "touch ${args.split('--arrow')[1].trim().split(' ')[0]}" : ""

    """
    ${arrow_file}
    touch ${prefix}.txt
    """
}
