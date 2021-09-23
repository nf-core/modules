#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BEDTOOLS_GENOMECOVSCALE } from '../../../../modules/bedtools/genomecovscale/main.nf' addParams( options: [suffix: '_out'] )

workflow test_bedtools_genomecovscale {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true),
        1
    ]

    sizes = file('dummy_chromosome_sizes')

    BEDTOOLS_GENOMECOVSCALE ( input, sizes )
}

workflow test_bedtools_genomecovscale_nonbam {
    input = [
        [ id:'test'],
        file(params.test_data['sarscov2']['genome']['baits_bed'], checkIfExists: true),
        1
    ]

    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    BEDTOOLS_GENOMECOVSCALE ( input, sizes )
}
