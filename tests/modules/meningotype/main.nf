#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MENINGOTYPE } from '../../../modules/meningotype/main.nf'

workflow test_meningotype {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    MENINGOTYPE ( input )
}
