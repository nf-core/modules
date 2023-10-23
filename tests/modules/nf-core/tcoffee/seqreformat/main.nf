#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TCOFFEE_SEQREFORMAT as TCOFFEE_SEQREFORMAT_SIMID  } from '../../../../../modules/nf-core/tcoffee/seqreformat/main.nf'

workflow test_tcoffee_seqreformat {
    
    input = [
        [ id:'test' ], // meta map
        file("https://raw.githubusercontent.com/nf-core/test-datasets/multiplesequencealign/testdata/setoxin-ref.fa", checkIfExists: true)
    ]

    TCOFFEE_SEQREFORMAT_SIMID ( input )
}
