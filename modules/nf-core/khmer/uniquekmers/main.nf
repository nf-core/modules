// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process KHMER_UNIQUEKMERS {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/36678b9530d1bf8c9541b42ca2f6dc1a7cf18b06257e876ab47a0afe9a63a97c/data'
        : 'community.wave.seqera.io/library/khmer:3.0.0a3--ad96a7538919c2bb'}"

    input:
    tuple val(meta), path(fasta)
    val kmer_size

    output:
    tuple val(meta), path("*.report.txt")  , emit: report
    tuple val(meta), path("*.kmers.txt")   , emit: kmers
    tuple val("${task.process}"), val('khmer'), eval('unique-kmers.py --version 2>&1 | grep ^khmer | sed "s/^khmer //;s/ .*$//"'), emit: versions_khmer, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    unique-kmers.py \\
        -k $kmer_size \\
        -R ${prefix}.report.txt \\
        $args \\
        $fasta

    grep ^number ${prefix}.report.txt | sed 's/^.*:.[[:blank:]]//g' > ${prefix}.kmers.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.report.txt
    touch ${prefix}.kmers.txt
    """
}
