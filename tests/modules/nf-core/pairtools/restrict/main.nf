#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_RESTRICT } from '../../../../../modules/nf-core/pairtools/restrict/main.nf'

workflow test_pairtools_restrict {

    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['generic']['pairtools']['mock_4flip_pairs'], checkIfExists: true) ]
    frag = file(params.test_data['generic']['pairtools']['frag_bed'], checkIfExists: true)

    PAIRTOOLS_RESTRICT ( input, frag )
}
