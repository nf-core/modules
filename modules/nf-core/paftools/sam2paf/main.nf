process PAFTOOLS_SAM2PAF {
    tag "$meta.id"
    label 'process_low'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37671219cfd244eb9b33db9345d3543ffd83037419a1c57f4648aace493ec2c2/data' :
        'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), file("*.paf"), emit: paf
    tuple val("${task.process}"), val('paftools'), eval("paftools.js version"), topic: versions, emit: versions_paftools
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    samtools view -h ${bam} | paftools.js sam2paf - > ${prefix}.paf
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.paf
    """
}
