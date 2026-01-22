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
    def coverage_stats        = "${prefix}_debug/coverage_stats.json"
    def sequences_chrM        = "${prefix}_debug/sequences_chrM.fa"
    def custom_alignments_bam = "${prefix}_debug/mito_igv_custom/custom_alignments.bam"
    def custom_alignments_bai = "${prefix}_debug/mito_igv_custom/custom_alignments.bam.bai"
    def igv_session           = "${prefix}_debug/mito_igv_custom/custom_igv_session.xml"
    def custom_ref_fa         = "${prefix}_debug/mito_igv_custom/custom_reference.fa"
    def custom_ref_fai        = "${prefix}_debug/mito_igv_custom/custom_reference.fa.fai"
    def custom_regions        = "${prefix}_debug/mito_igv_custom/custom_regions.bed"

    if (include_debug_output) {
        make_dirs = "mkdir -p ${prefix}_debug/mito_igv_custom"
        touch_debug_files = "touch ${coverage_stats} ${sequences_chrM} ${custom_alignments_bam} ${custom_alignments_bai} ${igv_session} ${custom_ref_fa} ${custom_ref_fai} ${custom_regions}"
    }
    else {
        make_dirs = ''
        touch_debug_files = ''
    }

    """
    echo $args

    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${touch_hap_stats}
    ${make_dirs}
    ${touch_debug_files}
    """
}
