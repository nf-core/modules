#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_INPUT1 } from '../../../../modules/gatk4/collectreadcounts/main.nf'
include { GATK4_COLLECTREADCOUNTS as GATK4_COLLECTREADCOUNTS_INPUT2 } from '../../../../modules/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../../modules/gatk4/determinegermlinecontigploidy/main.nf'

workflow test_gatk4_determinegermlinecontigploidy_cohort {
    
    input1 = [
        [ id:'test1' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
    ]
    input2 = [
        [ id:'test2' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true),
    ]
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    dict = file(params.test_data['homo_sapiens']['genome']['genome_dict'], checkIfExists: true)
    model = []
    priors =  file(params.test_data['homo_sapiens']['genome']['genome_ploidy_priors_tsv'], checkIfExists: true)
    
    GATK4_COLLECTREADCOUNTS_INPUT1 ( input1, fasta, fai, dict )
    GATK4_COLLECTREADCOUNTS_INPUT2 ( input2, fasta, fai, dict )
    input = GATK4_COLLECTREADCOUNTS_INPUT1.out.tsv.collect{ it[1] }
            .combine(GATK4_COLLECTREADCOUNTS_INPUT2.out.tsv.collect{ it[1] })
            .map { file1, file2 -> return [[id:"test"], [file1, file2]]}
    GATK4_DETERMINEGERMLINECONTIGPLOIDY ( input, model, priors )
}
