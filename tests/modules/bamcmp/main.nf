#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAMCMP } from '../../../modules/bamcmp/main.nf' addParams( options: [:] )

workflow test_bamcmp {
    
    input = [ [ id:'test'],
              file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true),
              file(params.test_data['sarscov2']['genome']['test2_bed'], checkIfExists: true)
            ]

    BAMCMP ( input )
}
