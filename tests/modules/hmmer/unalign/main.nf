#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_UNALIGN } from '../../../../modules/hmmer/unalign/main.nf'

workflow test_hmmer_unalign {
    
    input = [
        [ id:'test' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/hmmer/PF14720_seed.alnfaa.gz', checkIfExists: true)
    ]

    HMMER_UNALIGN ( input )
}
