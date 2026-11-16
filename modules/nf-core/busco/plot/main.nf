// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BUSCO_PLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd93359cce52c80ad33d58e18a694f37dc8aee9ccf2530cce8c42f9ade5c1f3f/data'
        : 'community.wave.seqera.io/library/busco:6.1.0--c502c35c2ab05540' }"

    input:
    path short_summary_json, stageAs: 'busco/*'

    output:
    path '*.png', emit: png
    tuple val("${task.process}"), val('busco'), eval("busco --version | sed 's/BUSCO //'"), emit: versions_busco, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'busco_figure'
    """
    busco \\
        --plot busco \\
        ${args}

    mv ./busco/busco_figure.png ${prefix}.png
    """

    stub:
    def prefix = task.ext.prefix ?: 'busco_figure'
    """
    touch ${prefix}.png
    """
}
