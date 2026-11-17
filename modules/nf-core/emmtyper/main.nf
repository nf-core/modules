// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process EMMTYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0ff2de4d7146b9da64cc2100c8fa40c0fbf52327affe11b0a45f17e395d2d44d/data'
        : 'community.wave.seqera.io/library/emmtyper:0.2.0--cc7e9aa35421a631'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('emmtyper'), eval('python -c "import emmtyper; print(emmtyper.__version__)"'), topic: versions, emit: versions_emmtyper

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    emmtyper \\
        $args \\
        $fasta \\
        > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
