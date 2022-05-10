#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GAMMA } from '../../../modules/gamma/main.nf'

workflow test_gamma {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    ]

    db = [ file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true) ]

    GAMMA ( input, db )
}
