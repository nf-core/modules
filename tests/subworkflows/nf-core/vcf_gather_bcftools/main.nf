#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VCF_GATHER_BCFTOOLS } from '../../../../subworkflows/nf-core/vcf_gather_bcftools/main.nf'

workflow test_vcf_gather_bcftools {

    input = Channel.of([
        [id:'test_1', sample:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ],
    [
        [id:'test_2', sample:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:'test_1', sample:'test'],
        2
    ],
    [
        [id:'test_2', sample:'test'],
        2
    ])


    VCF_GATHER_BCFTOOLS (
        input,
        scatter,
        [],
        'sample',
        true
    )
}

workflow test_vcf_gather_bcftools_no_meta {

    input = Channel.of([
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_ann_vcf_gz_tbi'], checkIfExists: true)
    ],
    [
        [id:'test'],
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_haplotc_vcf_gz_tbi'], checkIfExists: true)
    ])

    scatter = Channel.of([
        [id:'test'],
        2
    ],
    [
        [id:'test'],
        2
    ])

    bed = file(params.test_data['homo_sapiens']['genome']['genome_bed'], checkIfExists: true)


    VCF_GATHER_BCFTOOLS (
        input,
        scatter,
        bed,
        [],
        true
    )
}
