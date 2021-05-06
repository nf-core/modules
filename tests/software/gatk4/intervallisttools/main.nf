#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_INTERVALLISTTOOLS } from '../../../../software/gatk4/intervallisttools/main.nf' addParams( options: [:] )

workflow test_gatk4_intervallisttools {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]

    GATK4_INTERVALLISTTOOLS ( input )
}
