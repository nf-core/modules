#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AGRVATE } from '../../../modules/agrvate/main.nf' addParams( options: ["args": "--mummer"] )

workflow test_agrvate {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    AGRVATE ( input )
}
