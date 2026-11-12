// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BWAMEM2_INDEX {
    tag "$fasta"
    // NOTE Requires 28N GB memory where N is the size of the reference sequence, floor of 280M
    // source: https://github.com/bwa-mem2/bwa-mem2/issues/9
    // For bgzipped inputs this is scaled by a factor of 4.
    memory { 280.MB * Math.ceil(fasta.size() / 10000000) * task.attempt * (fasta.endsWith("gz") ? 4 : 1) }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e05ce34b46ad42810eb29f74e4e304c0cb592b2ca15572929ed8bbaee58faf01/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    tuple val(meta), path(fasta)
    val index_dir_name

    output:
    tuple val(meta), path("${outdir}"), emit: index
    tuple val("${task.process}"), val('bwamem2'), eval('bwa-mem2 version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_bwamem2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outdir = index_dir_name ?: 'bwamem2'
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    """
    mkdir ${outdir}
    bwa-mem2 \\
        index \\
        $args \\
        -p ${outdir}/${prefix} \\
        $fasta
    """

    stub:
    outdir = index_dir_name ?: 'bwamem2'
    def prefix = task.ext.prefix ?: "${fasta}"

    """
    mkdir ${outdir}
    touch ${outdir}/${prefix}.0123
    touch ${outdir}/${prefix}.ann
    touch ${outdir}/${prefix}.pac
    touch ${outdir}/${prefix}.amb
    touch ${outdir}/${prefix}.bwt.2bit.64
    """
}
