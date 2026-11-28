// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ARTIC_GUPPYPLEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f6019bce7e1652ddff248bbfe2a6d9fc1c668137b731cd5ed772e690a7d5ddca/data' :
        'community.wave.seqera.io/library/artic:1.12.0--ac55155577657f92' }"

    input:
    tuple val(meta), path(fastq_dir)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('artic'), eval('artic -v 2>&1 | sed "s/^.*artic //; s/ .*$//"')   , emit: versions_artic, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic \\
        guppyplex \\
        ${args} \\
        --directory ${fastq_dir} \\
        --output ${prefix}.fastq

    pigz -p ${task.cpus} *.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip() { pigz; }
    echo "" | gzip > ${prefix}.fastq.gz
    """
}
