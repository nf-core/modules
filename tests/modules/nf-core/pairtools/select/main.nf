#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_SELECT } from '../../../../../modules/nf-core/pairtools/select/main.nf'

workflow test_pairtools_select {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['generic']['pairtools']['mock_pairsam'], checkIfExists: true) ]

    PAIRTOOLS_SELECT ( input )
}
