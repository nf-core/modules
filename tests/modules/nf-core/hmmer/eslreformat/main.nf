#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMALIGN                                 } from '../../../../../modules/nf-core/hmmer/hmmalign/main.nf'
include { HMMER_ESLREFORMAT as HMMER_ESLREFORMAT_AFA     } from '../../../../../modules/nf-core/hmmer/eslreformat/main.nf'
include { HMMER_ESLREFORMAT as HMMER_ESLREFORMAT_PHYLIP  } from '../../../../../modules/nf-core/hmmer/eslreformat/main.nf'
include { HMMER_ESLREFORMAT as HMMER_ESLREFORMAT_UNALIGN } from '../../../../../modules/nf-core/hmmer/eslreformat/main.nf'

workflow test_hmmer_eslreformat_afa {
    
    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/e_coli_k12_16s.fna.gz')      // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/bac.16S_rRNA.hmm.gz')

    HMMER_HMMALIGN ( input, hmm )

    HMMER_ESLREFORMAT_AFA ( HMMER_HMMALIGN.out.sthlm )
}

workflow test_hmmer_eslreformat_phylip {
    
    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/e_coli_k12_16s.fna.gz')      // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/bac.16S_rRNA.hmm.gz')

    HMMER_HMMALIGN ( input, hmm )

    HMMER_ESLREFORMAT_PHYLIP ( HMMER_HMMALIGN.out.sthlm )
}

workflow test_hmmer_eslreformat_unalign {
    
    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/e_coli_k12_16s.fna.gz')      // Change to params.test_data syntax after the data is included in ./tests/config/test_data.config
    ]

    hmm   = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/delete_me/hmmer/bac.16S_rRNA.hmm.gz')

    HMMER_HMMALIGN ( input, hmm )

    HMMER_ESLREFORMAT_UNALIGN ( HMMER_HMMALIGN.out.sthlm )
}
