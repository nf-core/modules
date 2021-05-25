#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BIGWIGAVERAGEOVERBED } from '../../../../software/ucsc/bigwigaverageoverbed/main.nf' addParams( options: [:] )

workflow test_ucsc_bigwigaverageoverbed {
    input = [ [ id: 'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true ) ]
            ]
    bigwig = file(params.test_data['sarscov2']['illumina']['test_bigwig'], checkIfExists: true)

    UCSC_BIGWIGAVERAGEOVERBED ( input, bigwig )
}
