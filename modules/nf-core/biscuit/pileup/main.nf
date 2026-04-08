process BISCUIT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/2790659e7b390a29fa72483e7eed0101298dd2959dce2f8f3744b4f93517d323/data':
        'community.wave.seqera.io/library/biscuit_htslib:b5833b773bc32c9c' }"

    input:
    tuple val(meta), path(normal_bams), path(normal_bais), path(tumor_bam), path(tumor_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('biscuit'), eval("biscuit version |& sed '1!d; s/^.*BISCUIT Version: //'"), emit: versions_biscuit, topic: versions
    tuple val("${task.process}"), val("bgzip"), eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"')    , emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biscuit_cpus = [(task.cpus * 0.9) as int, 1].max()
    def bgzip_cpus = (task.cpus - biscuit_cpus < 1) ? biscuit_cpus : (task.cpus - biscuit_cpus)
    if ( tumor_bam != [] && normal_bams.toList().size() > 1 ) error "[BISCUIT_PILEUP] error: Tumor BAM provided with more than one normal BAM"
    if ( tumor_bam.toList().size() > 1 ) error "[BISCUIT_PILEUP] error: more than one tumor BAM provided"
    input = ( tumor_bam==[] ) ? "${normal_bams}" : "-S -T ${tumor_bam} -I ${normal_bams}"
    """
    ln -sf \$(readlink $fasta) $index/$fasta

    biscuit pileup \\
        -@ $biscuit_cpus \\
        $args \\
        $index/$fasta \\
        $input \\
        | bgzip -@ $bgzip_cpus $args2 > ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    """
}
