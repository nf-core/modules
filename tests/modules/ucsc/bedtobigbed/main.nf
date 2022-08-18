#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { UCSC_BEDTOBIGBED } from '../../../../modules/ucsc/bedtobigbed/main.nf'

workflow test_ucsc_bedtobigbed {
    input = [ [ id: 'test' ], // meta map 
              [ file(params.test_data['sarscov2']['genome']['test_bed12'], checkIfExists: true ) ] 
            ]
    sizes = file(params.test_data['sarscov2']['genome']['genome_sizes'], checkIfExists: true)

    UCSC_BEDTOBIGBED ( input, sizes )
}
