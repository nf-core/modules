#!/usr/bin/env nextflow



include { NGMASTER } from '../../../modules/ngmaster/main.nf'

workflow test_ngmaster {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    NGMASTER ( input )
}
