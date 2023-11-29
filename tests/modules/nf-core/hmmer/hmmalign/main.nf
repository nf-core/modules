#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMALIGN as HMMER_HMMALIGN        } from '../../../../../modules/nf-core/hmmer/hmmalign/main.nf'
include { HMMER_HMMALIGN as HMMER_HMMALIGN_PREFIX } from '../../../../../modules/nf-core/hmmer/hmmalign/main.nf'

workflow test_hmmer_hmmalign {

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/e_coli_k12_16s.fna.gz')      // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/bac.16S_rRNA.hmm.gz')

    HMMER_HMMALIGN ( input, hmm )
}

workflow test_hmmer_hmmalign_prefix {

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/e_coli_k12_16s.fna.gz')      // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/bac.16S_rRNA.hmm.gz')

    HMMER_HMMALIGN_PREFIX ( input, hmm )
}
