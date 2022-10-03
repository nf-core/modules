#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PASTY } from '../../../../modules/nf-core/pasty/main.nf'

workflow test_pasty {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    PASTY ( input )
}
