process MITORSAW_HAPLOTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitorsaw:0.2.7--h9ee0642_0':
        'biocontainers/mitorsaw:0.2.7--h9ee0642_0' }"

    input:

    tuple val(meta), path(bam), path(bai)    // channel: [ val(meta), path(bam), path(bai) ]
    tuple val(meta2), path(fasta), path(fai) // channel: [ val(meta2), path(fasta), path(fai) ]
    val(include_hap_stats)                   // boolean:   [ true | false ], optional: true
    val(include_debug_output)                // boolean: [ true | false ], optional: true

    output:

    tuple val(meta), path("*${prefix}.vcf.gz"), path("*${prefix}.vcf.gz.tbi"),                                                                         emit: vcf
    tuple val(meta), path("${prefix}.json"),                                                                                                           emit: stats,             optional: true
    tuple val(meta), path("${prefix}_debug/coverage_stats.json"),                                                                                      emit: coverage_stats,    optional: true
    tuple val(meta), path("${prefix}_debug/sequences_chrM.fa"),                                                                                        emit: sequences_chrM,    optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_alignments.bam"), path("${prefix}_debug/mito_igv_custom/custom_alignments.bam.bai"), emit: custom_alignments, optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_igv_session.xml"),                                                                   emit: igv_session,       optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_reference.fa"), path("${prefix}_debug/mito_igv_custom/custom_reference.fa.fai"),     emit: custom_ref,        optional: true
    tuple val(meta), path("${prefix}_debug/mito_igv_custom/custom_regions.bed"),                                                                       emit: custom_regions,    optional: true

    tuple val("${task.process}"), val('mitorsaw'), eval("mitorsaw --version | cut -d' ' -f2 | cut -d'-' -f1"), topic: versions, emit: versions_mitorsaw

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
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
