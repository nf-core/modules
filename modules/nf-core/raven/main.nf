// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process RAVEN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ec/ece0d4ba37373e4fcf84fd856b4e7f3b84ecc4394dabeed7928a0c9a9b3b92bb/data'
        : 'community.wave.seqera.io/library/raven-assembler:1.6.1--2388e7e21e19cd95'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    tuple val("${task.process}"), val('raven'), eval('raven --version'), emit: versions_raven, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # run tool
    raven \\
        -t $task.cpus \\
        --graphical-fragment-assembly ${prefix}.gfa \\
        $args \\
        $reads | \\
        gzip -c > ${prefix}.fasta.gz

    # compress assembly graph
    gzip -c ${prefix}.gfa > ${prefix}.gfa.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fasta.gz
    echo "" | gzip > ${prefix}.gfa.gz
    """
}
