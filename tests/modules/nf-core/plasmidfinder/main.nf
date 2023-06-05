#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PLASMIDFINDER } from '../../../../modules/nf-core/plasmidfinder/main.nf'

workflow test_plasmidfinder {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]

    PLASMIDFINDER ( input )
}
