#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PARABRICKS_HAPLOTYPECALLER } from '../../../../../modules/nf-core/parabricks/haplotypecaller/main.nf'
include { PARABRICKS_HAPLOTYPECALLER as PARABRICKS_HAPLOTYPECALLER_GVCF } from '../../../../../modules/nf-core/parabricks/haplotypecaller/main.nf'

workflow test_parabricks_haplotypecaller {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        [],
        []
    ]
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PARABRICKS_HAPLOTYPECALLER ( input, fasta )
}

workflow test_parabricks_haplotypecaller_intervals {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)      
    ]

    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PARABRICKS_HAPLOTYPECALLER ( input, fasta )
}

workflow test_parabricks_haplotypecaller_gvcf {

    input = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_recalibrated_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)
    ]
    
    fasta = [
        [ id:'test'],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    PARABRICKS_HAPLOTYPECALLER_GVCF ( input, fasta )
}
