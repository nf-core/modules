// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWA_INDEX {
    tag "$fasta"
    // NOTE requires 5.37N memory where N is the size of the database
    // source: https://bio-bwa.sourceforge.net/bwa.shtml#8
    memory { 7.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(fasta)
    val input_dir

    output:
    tuple val(meta), path("${outdir}"), emit: index
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa

    when:
    task.ext.when == null || task.ext.when

    script:
    outdir = input_dir ?: 'bwa'
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir ${outdir}
    bwa \\
        index \\
        $args \\
        -p ${outdir}/${prefix} \\
        $fasta
    """

    stub:
    outdir = input_dir ?: 'bwa'
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir ${outdir}
    touch ${outdir}/${prefix}.amb
    touch ${outdir}/${prefix}.ann
    touch ${outdir}/${prefix}.bwt
    touch ${outdir}/${prefix}.pac
    touch ${outdir}/${prefix}.sa
    """
}
