process HIPHASE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d418cdbaff565e9c563c441c72d480f2605bb529712dd026068f1d0a7b246617/data'
        : 'community.wave.seqera.io/library/hiphase:1.5.0--f36e5874e9287052'}"

    input:
    tuple val(meta), path(bams), path(bais), path(snv), path(snv_index), path(sv), path(sv_index), val(samples)
    tuple val(meta2), path(fasta), path(fai)
    val output_bam

    output:
    tuple val(meta), path("*_snv_phased.vcf.gz"), emit: vcfs, optional: true
    tuple val(meta), path("*_sv_phased.vcf.gz"), emit: sv_vcfs, optional: true
    tuple val(meta), path("*_snv_phased.vcf.gz.{tbi,csi}"), emit: vcfs_indexes, optional: true
    tuple val(meta), path("*_sv_phased.vcf.gz.{tbi,csi}"), emit: sv_vcfs_indexes, optional: true
    tuple val(meta), path("*.summary.tsv"), emit: summary_tsv, optional: true
    tuple val(meta), path("*.summary.csv"), emit: summary_csv, optional: true
    tuple val(meta), path("*.blocks.tsv"), emit: blocks_tsv, optional: true
    tuple val(meta), path("*.blocks.csv"), emit: blocks_csv, optional: true
    tuple val(meta), path("*.stats.tsv"), emit: stats_tsv, optional: true
    tuple val(meta), path("*.stats.csv"), emit: stats_csv, optional: true
    tuple val(meta), path("*.haplotag.tsv"), emit: haplotag_tsv, optional: true
    tuple val(meta), path("*.haplotag.csv"), emit: haplotag_csv, optional: true
    tuple val(meta), path("*.bam"), emit: bams, optional: true
    tuple val(meta), path("*.bam.{bai,csi}"), emit: bams_indexes, optional: true
    tuple val("${task.process}"), val('hiphase'), eval("hiphase --version |& sed '1!d ; s/hiphase //'"), emit: versions_hiphase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // This is set up to allow for co-phasing SNVs and SVs, each input with their own VCF and output separately.
    // The VCF files can contain multiple samples, with multiple BAM files, and the samples to be phased can be specified with the sample argument.
    def snv_args = snv ? "--vcf ${snv} --output-vcf ${prefix}_snv_phased.vcf.gz" : ''
    def sv_args = sv ? "--vcf ${sv} --output-vcf ${prefix}_sv_phased.vcf.gz" : ''
    // Cannot use prefix for bam outputs as we can have multiple output BAM files leading to name collisions.
    def bam_args = bams
        .collectMany { file ->
            ["--bam", file, output_bam ? '--output-bam' : '', output_bam ? "${file.baseName}_haplotagged.bam" : '']
        }
        .join(" ")

    def sample_args = samples ? samples.collect { sample -> "--sample-name ${sample}" }.join(" ") : ''

    """
    hiphase \
        ${args} \
        --threads ${task.cpus} \\
        --reference ${fasta} \\
        ${sample_args} \\
        ${bam_args} \\
        ${snv_args} \\
        ${sv_args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    def bam_index_format = args.contains('--csi-index') ? 'csi' : 'bai'
    def vcf_index_format = args.contains('--csi-index') ? 'csi' : 'tbi'
    def touch_bams = output_bam ? bams.collect { file -> "touch ${prefix}_${file.baseName}_haplotagged.bam" }.join("\n") : ''
    def touch_bams_indexes = output_bam ? bams.collect { file -> "touch ${prefix}_${file.baseName}_haplotagged.bam.${bam_index_format}" }.join("\n") : ''
    def touch_snv_vcf = snv ? "echo | gzip > ${prefix}_snv_phased.vcf.gz" : ''
    def touch_snv_vcf_index = snv ? "touch ${prefix}_snv_phased.vcf.gz.${vcf_index_format}" : ''
    def touch_sv_vcf = sv ? "echo | gzip > ${prefix}_sv_phased.vcf.gz" : ''
    def touch_sv_vcf_index = sv ? "touch ${prefix}_sv_phased.vcf.gz.${vcf_index_format}" :''
    """
    ${touch_bams}
    ${touch_bams_indexes}
    ${touch_snv_vcf}
    ${touch_snv_vcf_index}
    ${touch_sv_vcf}
    ${touch_sv_vcf_index}
    """
}
