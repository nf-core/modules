#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_QUERY } from '../../../../modules/bcftools/query/main.nf'

workflow test_bcftools_query {

    regions = []
    targets = []
    samples = []

    input = [ [ id:'out' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    BCFTOOLS_QUERY ( input, regions, targets, samples )
}

workflow test_bcftools_query_with_optional_files {

    regions = file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true)
    targets = file(params.test_data['sarscov2']['illumina']['test2_vcf_targets_tsv_gz'], checkIfExists: true)
    samples = []

    input = [ [ id:'out' ], // meta map
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
             file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]

    BCFTOOLS_QUERY ( input, regions, targets, samples )
}
