// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NGSCHECKMATE_VAFNCM {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/643c1bd726dd16bf3a634ed3ffc0bfa573b5b254da7e0450a1e910bd1f4619d6/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_ngscheckmate:f79d3ada00495af3'}"

    input:
    tuple val(meta), path(vafs)

    output:
    tuple val(meta), path("*.pdf")             , emit: pdf, optional: true
    tuple val(meta), path("*_corr_matrix.txt") , emit: corr_matrix
    tuple val(meta), path("*_all.txt")         , emit: all
    tuple val(meta), path("*_matched.txt")     , emit: matched
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # tool has a bug where it misses the final file, so add a dummy one.
    cp ${vafs[0]} zzzzzz.vaf
    vaf_ncm.py -I . -O . -N ${prefix} $args

    # remove the existence of the dummy file
    rm zzzzzz.vaf
    sed -i.bak "/zzzzzz/d" ${prefix}_all.txt

    # generate a file with all the samples that do match, for consistency with the bam mode (ngscheckmate/ncm)
    sed "/unmatched/d" ${prefix}_all.txt > ${prefix}_matched.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf
    """
}
