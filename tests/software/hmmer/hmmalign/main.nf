#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HMMER_HMMALIGN } from '../../../../software/hmmer/hmmalign/main.nf' addParams( options: [:] )

workflow test_hmmer_hmmalign {
    
//    input = [ [ id:'test' ], // meta map
//              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    input = [
        [ id:'test' ], // meta map
        file('https://raw.githubusercontent.com/erikrikarddaniel/test-datasets/modules/data/delete_me/e_coli_k12_16s.fna')
    ]

    hmm   = file('https://raw.githubusercontent.com/erikrikarddaniel/test-datasets/modules/data/delete_me/bac.16S_rRNA.hmm')

    HMMER_HMMALIGN ( input, hmm )
}

// vim:sw=4
