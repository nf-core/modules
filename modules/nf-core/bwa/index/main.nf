process BWA_INDEX {
    tag "$fasta"
    // NOTE requires 5.37N memory where N is the size of the database
    // source: https://bio-bwa.sourceforge.net/bwa.shtml#8
    memory { 6.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/"), emit: index
    tuple val("${task.process}"), val('bwa'), eval('bwa 2>&1 | sed -n "s/^Version: //p"'), topic: versions, emit: versions_bwa

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    def genome = "${fasta.baseName}"
    """
    mkdir ${prefix}
    bwa \\
        index \\
        $args \\
        -p ${prefix}/${genome} \\
        $fasta
    """

    stub:
    prefix = task.ext.prefix ?: meta.id ?: "${fasta.baseName}"
    def genome = "${fasta.baseName}"
    """
    mkdir ${prefix}
    touch ${prefix}/${genome}.amb
    touch ${prefix}/${genome}.ann
    touch ${prefix}/${genome}.bwt
    touch ${prefix}/${genome}.pac
    touch ${prefix}/${genome}.sa
    """
}
