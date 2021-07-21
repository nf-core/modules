#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../modules/bedtools/genomecov/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_genomecov {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    sizes = file('dummy_chromosome_sizes')
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}

workflow test_bedtools_genomecov_nonbam {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)
    ]

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)
    extension = 'txt'

    BEDTOOLS_GENOMECOV ( input, sizes, extension )
}
