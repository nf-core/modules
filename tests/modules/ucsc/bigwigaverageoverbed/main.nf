#!/usr/bin/env nextflow



include { UCSC_BIGWIGAVERAGEOVERBED } from '../../../../modules/ucsc/bigwigaverageoverbed/main.nf'

workflow test_ucsc_bigwigaverageoverbed {
    input = [
        [ id: 'test' ], // meta map
        [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true ) ]
    ]
    bigwig = file(params.test_data['sarscov2']['illumina']['test_bigwig'], checkIfExists: true)

    UCSC_BIGWIGAVERAGEOVERBED ( input, bigwig )
}
