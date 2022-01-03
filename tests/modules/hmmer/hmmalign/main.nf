#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMALIGN } from '../../../../modules/hmmer/hmmalign/main.nf'

workflow test_hmmer_hmmalign {

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/erikrikarddaniel/test-datasets/modules/data/delete_me/e_coli_k12_16s.fna')      // Change to params.test_data syntax after the data is included in tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/erikrikarddaniel/test-datasets/modules/data/delete_me/bac.16S_rRNA.hmm')

    HMMER_HMMALIGN ( input, hmm )
}
