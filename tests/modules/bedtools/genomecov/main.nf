#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../modules/bedtools/genomecov/main.nf'

workflow test_bedtools_genomecov_noscale {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        1
    ]

    sizes = []
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}

workflow test_bedtools_genomecov_nonbam_noscale {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true),
        1
    ]

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}

workflow test_bedtools_genomecov_scale {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        0.5
    ]

    sizes = file('dummy_chromosome_sizes')
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}

workflow test_bedtools_genomecov_nonbam_scale {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true),
        0.5
    ]

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}
