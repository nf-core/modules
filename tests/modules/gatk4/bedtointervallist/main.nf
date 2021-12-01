#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BEDTOINTERVALLIST } from '../../../../modules/gatk4/bedtointervallist/main.nf'

workflow test_gatk4_bedtointervallist {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) ]
            ]
    dict = file(params.test_data['sarscov2']['genome']['genome_dict'], checkIfExists: true)

    GATK4_BEDTOINTERVALLIST ( input, dict )
}
