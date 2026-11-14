// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process GRAPHMAP2_INDEX {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2c/2cc3104dc1685c517113b64d39a936d04977410994ec5982ef22fd2b76199d8f/data'
        : 'community.wave.seqera.io/library/graphmap:0.6.4--4c694751d262fa1f'}"

    input:
    path fasta

    output:
    path "*.gmidx", emit: index
    tuple val("${task.process}"), val('graphmap2'), eval("graphmap2 2>&1 | sed -n 's/Version: v//p'"), emit: versions_graphmap2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    graphmap2 \\
        align \\
        -t ${task.cpus} \\
        -I \\
        ${args} \\
        -r ${fasta}
    """

    stub:

    """
    touch ${fasta}.gmidx
    """
}
