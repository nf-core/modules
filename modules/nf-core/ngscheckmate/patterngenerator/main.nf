// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NGSCHECKMATE_PATTERNGENERATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/643c1bd726dd16bf3a634ed3ffc0bfa573b5b254da7e0450a1e910bd1f4619d6/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_ngscheckmate:f79d3ada00495af3'}"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bowtie_index)

    output:
    tuple val(meta), path("*.pt"), emit: pt
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')

    makesnvpattern.pl ${bed} ${fasta} \$INDEX . ${prefix}
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$fasta" == "${prefix}.fasta") error "makesnvpattern.pl generates a fasta file with the same name as the input fasta, use \"task.ext.prefix\" to disambiguate!"

    """
    touch ${prefix}.pt
    """
}
