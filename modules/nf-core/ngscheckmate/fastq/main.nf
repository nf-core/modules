// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NGSCHECKMATE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/643c1bd726dd16bf3a634ed3ffc0bfa573b5b254da7e0450a1e910bd1f4619d6/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_ngscheckmate:f79d3ada00495af3'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(snp_pt)

    output:
    tuple val(meta), path("*.vaf"), emit: vaf
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq2command  = ( reads instanceof List && reads.size() == 2 ) ? " -2 ${reads[1]} " : ""

    """
    ngscheckmate_fastq -1 ${reads[0]} $fastq2command ${snp_pt} -p ${task.cpus} $args > ${prefix}.vaf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vaf
    """
}
