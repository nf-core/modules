process MITORSAW_HAPLOTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitorsaw:0.2.7--h9ee0642_0':
        'biocontainers/mitorsaw:0.2.7--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    val(include_hap_stats)
    val(include_debug_output)

    output:
    tuple val(meta), path("*${prefix}.vcf.gz"),                                         emit: vcf
    tuple val(meta), path("*${prefix}.vcf.gz.tbi"),                                     emit: tbi
    tuple val(meta), path("${prefix}.json"),                                            emit: stats,                 optional: true
    tuple val(meta), path("${prefix}_debug/coverage_stats.json"),                       emit: coverage_stats,        optional: true
    tuple val(meta), path("${prefix}_debug/sequences_chrM.fa"),                         emit: sequences_chrM,        optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_alignments.bam"),     emit: custom_alignments_bam, optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_alignments.bam.bai"), emit: custom_alignments_bai, optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_igv_session.xml"),    emit: igv_session,           optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_reference.fa"),       emit: custom_ref_fa,         optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_reference.fa.fai"),   emit: custom_ref_fai,        optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_regions.bed"),        emit: custom_regions,        optional: true

    tuple val("${task.process}"), val('mitorsaw'), eval("mitorsaw --version | sed 's/mitorsaw //'"), topic: versions, emit: versions_mitorsaw

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def output_hap_stats = include_hap_stats    ? "--output-hap-stats ${prefix}.json" : ''
    def debug_output     = include_debug_output ? "--output-debug ${prefix}_debug"    : ''

    """
    mitorsaw haplotype \\
        $args \\
        --reference ${fasta} \\
        --bam ${bam} \\
        --output-vcf ${prefix}.vcf.gz \\
        ${output_hap_stats} \\
        ${debug_output}
    """

    stub:
    def args                  = task.ext.args        ?: ''
    prefix                    = task.ext.prefix      ?: "${meta.id}"
    def touch_hap_stats       = include_hap_stats    ? "touch ${prefix}.json" : ''
    def debug_dir             = "${prefix}_debug"
    def igv_dir               = "${debug_dir}/mito_igv_custom"
    def mkdir_debug           = include_debug_output ? "mkdir -p ${igv_dir}"  : ''
    def debug_files           = [
        "${debug_dir}/coverage_stats.json",
        "${debug_dir}/sequences_chrM.fa",
        "${igv_dir}/custom_alignments.bam",
        "${igv_dir}/custom_alignments.bam.bai",
        "${igv_dir}/custom_igv_session.xml",
        "${igv_dir}/custom_reference.fa",
        "${igv_dir}/custom_reference.fa.fai",
        "${igv_dir}/custom_regions.bed"
    ]
    def touch_debug_files = include_debug_output ? "touch ${debug_files.join(' ')}" : ''
    """
    echo $args

    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${touch_hap_stats}
    ${mkdir_debug}
    ${touch_debug_files}
    """
}
