// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process NGSCHECKMATE_NCM {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/643c1bd726dd16bf3a634ed3ffc0bfa573b5b254da7e0450a1e910bd1f4619d6/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_ngscheckmate:f79d3ada00495af3'}"

    input:
    tuple val(meta), path(files)
    tuple val(meta2), path(snp_bed)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path("*_corr_matrix.txt"), emit: corr_matrix
    tuple val(meta), path("*_matched.txt"),     emit: matched
    tuple val(meta), path("*_all.txt"),         emit: all
    tuple val(meta), path("*.pdf"),             emit: pdf, optional: true
    tuple val(meta), path("*.vcf"),             emit: vcf, optional: true
    tuple val("${task.process}"), val('ngscheckmate'), eval("ncm.py --help | sed '7!d;s/.* v//g'"), topic: versions, emit: versions_ngscheckmate
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unzip = files.any {file ->  file.toString().endsWith(".vcf.gz") }
    """
    if ${unzip}
    then
        for VCFGZ in *.vcf.gz; do
            gunzip -cdf \$VCFGZ > \$( basename \$VCFGZ .gz );
        done
    fi

    NCM_REF="./"${fasta} ncm.py -d . -bed ${snp_bed} -O . -N ${prefix} ${args}

    if ${unzip}
    then
        rm -f *.vcf  # clean up decompressed vcfs
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_output_corr_matrix.txt
    touch ${prefix}_matched.txt
    touch ${prefix}_all.txt
    touch ${prefix}.pdf
    """
}
