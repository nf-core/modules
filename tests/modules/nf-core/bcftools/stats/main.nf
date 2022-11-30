#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BCFTOOLS_STATS } from '../../../../../modules/nf-core/bcftools/stats/main.nf'

workflow test_bcftools_stats {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions = []
    targets = []
    samples = []

    BCFTOOLS_STATS ( input, regions, targets, samples )
}

workflow test_bcftools_stats_regions {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz_tbi'], checkIfExists: true)]
    regions = file(params.test_data['sarscov2']['illumina']['test3_vcf_gz'], checkIfExists: true)
    targets = []
    samples = []

    BCFTOOLS_STATS ( input, regions, targets, samples )
}

workflow test_bcftools_stats_targets {
    input = [ [ id:'test' ], // meta map
            file(params.test_data['sarscov2']['illumina']['test_vcf_gz'], checkIfExists: true),
            []]
    regions = []
    targets = file(params.test_data['sarscov2']['illumina']['test2_vcf_targets_tsv_gz'], checkIfExists: true)
    samples = []

    BCFTOOLS_STATS ( input, regions, targets, samples )
}
