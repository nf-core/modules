#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include { GATK4_PREPROCESSINTERVALS } from '../../../../modules/gatk4/preprocessintervals/main.nf'
// include { GATK4_ANNOTATEINTERVALS } from '../../../../modules/gatk4/annotateintervals/main.nf'
// include { GATK4_COLLECTREADCOUNTS } from '../../../../modules/gatk4/collectreadcounts/main.nf'
include { GATK4_FILTERINTERVALS } from '../../../../modules/gatk4/filterintervals/main.nf'

workflow test_gatk4_filterintervals {
    
    // bed = file(params.test_data['homo_sapiens']['genome']['genome_blacklist_interval_bed'], checkIfExists: true)
    // fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    // fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    // dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    // bam = [
    //     [ id:'test' ], // meta map
    //     file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
    //     file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
    // ]
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['genome']['genome_preprocessed_count_tsv'], checkIfExists: true)
    ]
    preprocessed_interval = file(params.test_data['homo_sapiens']['genome']['genome_preprocessed_interval_bed'], checkIfExists: true)
    annotated_interval = file(params.test_data['homo_sapiens']['genome']['genome_annotated_interval_tsv'], checkIfExists: true)

    // GATK4_PREPROCESSINTERVALS ( bed, fasta, fai, dict )
    // GATK4_ANNOTATEINTERVALS ( GATK4_PREPROCESSINTERVALS.out.interval_list, fasta, fai, dict  )
    // GATK4_COLLECTREADCOUNTS ( bam, GATK4_PREPROCESSINTERVALS.out.interval_list, fasta, fai, dict )
    // GATK4_FILTERINTERVALS ( GATK4_COLLECTREADCOUNTS.out.tsv, GATK4_PREPROCESSINTERVALS.out.interval_list, GATK4_ANNOTATEINTERVALS.out.tsv )
    GATK4_FILTERINTERVALS ( input, preprocessed_interval, annotated_interval )
}
