#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RTGTOOLS_ROCPLOT } from '../../../../../modules/nf-core/rtgtools/rocplot/main.nf'
include { RTGTOOLS_VCFEVAL } from '../../../../../modules/nf-core/rtgtools/vcfeval/main.nf'
include { UNTAR            } from '../../../../../modules/nf-core/untar/main.nf'

workflow test_rtgtools_rocplot {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)
    ]

    compressed_sdf = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    ]

    sdf = UNTAR( compressed_sdf ).untar

    RTGTOOLS_VCFEVAL ( input, sdf )

    RTGTOOLS_ROCPLOT ( RTGTOOLS_VCFEVAL.out.weighted_roc )
}


