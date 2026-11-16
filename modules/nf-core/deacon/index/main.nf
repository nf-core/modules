// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process DEACON_INDEX {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f628781900b19a962ce706a66d7b1d83c26acb5b08d5085947ae1bbea879b1f8/data'
        : 'community.wave.seqera.io/library/deacon:0.15.0--ac1f965b3a6ea82e'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.idx"), emit: index
    tuple val("${task.process}"), val('deacon'), eval('deacon --version | head -n1 | sed "s/deacon //g"'), emit: versions_deacon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    deacon \\
        index \\
        build \\
        --threads ${task.cpus} \\
        $args \\
        $fasta > ${prefix}.idx
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    touch ${prefix}.idx
    """
}
