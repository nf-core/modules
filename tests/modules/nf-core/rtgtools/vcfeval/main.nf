#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
moduleDir = launchDir

include { RTGTOOLS_VCFEVAL } from "$moduleDir/modules/nf-core/rtgtools/vcfeval/main.nf"
include { UNTAR } from "$moduleDir/modules/nf-core/untar/main.nf"

workflow test_rtgtools_vcfeval {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true),
    ]

    truth = [
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ]

    truth_regions = file(params.test_data['homo_sapiens']['genome']['genome_21_multi_interval_bed'], checkIfExists: true)

    evaluation_regions = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)

    compressed_sdf = [
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    ]

    sdf = UNTAR( compressed_sdf ).untar
        .map({
            meta, folder ->
                folder
        })


    RTGTOOLS_VCFEVAL ( input, truth, truth_regions, evaluation_regions, sdf )
}

workflow test_rtgtools_vcfeval_no_optional_inputs {

    input = [
        [ id:'test' ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        [],
    ]

    truth = [
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        []
    ]

    truth_regions = []

    evaluation_regions = []

    compressed_sdf = [
        [],
        file(params.test_data['homo_sapiens']['genome']['genome_21_sdf'])
    ]

    sdf = UNTAR( compressed_sdf ).untar
        .map({
            meta, folder ->
                [folder]
        })

    RTGTOOLS_VCFEVAL ( input, truth, truth_regions, evaluation_regions, sdf )
}
