#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMFETCH } from '../../../../../modules/nf-core/hmmer/hmmfetch/main.nf'

workflow test_hmmer_hmmfetch_key {
    
    input = file('https://raw.githubusercontent.com/tseemann/barrnap/master/db/arc.hmm')

    HMMER_HMMFETCH ( input, '16S_rRNA', [], [] )
}

workflow test_hmmer_hmmfetch_keyfile {
    
    input = file('https://raw.githubusercontent.com/tseemann/barrnap/master/db/arc.hmm')

    Channel
        .of('16S_rRNA', '23S_rRNA')
        .collectFile(name: 'keys.txt', newLine: true)
        .set { keyfile }

    HMMER_HMMFETCH ( input, [], keyfile, [] )
}

workflow test_hmmer_hmmfetch_index {
    
    input = file('https://raw.githubusercontent.com/tseemann/barrnap/master/db/arc.hmm')

    HMMER_HMMFETCH ( input, [], [], [] )
}
