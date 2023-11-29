#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_FLIP } from '../../../../../modules/nf-core/pairtools/flip/main.nf'

workflow test_pairtools_flip {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['generic']['pairtools']['mock_4flip_pairs'], checkIfExists: true) ]
    sizes = file(params.test_data['generic']['pairtools']['mock_chrom_sizes'], checkIfExists:true)

    PAIRTOOLS_FLIP ( input, sizes )
}
