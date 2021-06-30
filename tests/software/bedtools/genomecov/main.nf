#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOV } from '../../../../software/bedtools/genomecov/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_genomecov {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
            ]

    chromosome_sizes = file('dummy_chromosome_sizes')

    output_suffix = 'txt'

    BEDTOOLS_GENOMECOV ( input, chromosome_sizes, output_suffix )
}

workflow test_bedtools_genomecov_nonbam {
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true)
            ]

    chromosome_sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    output_suffix = 'txt'

    BEDTOOLS_GENOMECOV ( input, chromosome_sizes, output_suffix )
}