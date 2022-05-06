#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_VCFEVAL } from '../../../../modules/rtgtools/vcfeval/main.nf'

workflow test_rtgtools_vcfeval {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)        
    ]

    sdf = Channel.value(
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    )

    RTGTOOLS_VCFEVAL ( input, sdf )
}

workflow test_rtgtools_vcfeval_no_index {
    
    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        [],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)        
    ]

    sdf = Channel.value(
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    )

    RTGTOOLS_VCFEVAL ( input, sdf )
}
