#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BED12TOBIGBED } from '../../../../software/ucsc/bed12tobigbed/main.nf' addParams( options: [:] )

workflow test_ucsc_bed12tobigbed {
    input = [ [ id: 'test' ], // meta map 
              [ file(params.test_data['sarscov2']['genome']['test_bed12'], checkIfExists: true ) ] 
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BED12TOBIGBED ( input, sizes )
}
