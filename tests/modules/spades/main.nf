#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../../../software/spades/main.nf' addParams( spades_hmm: false ,options: ['args': '--rnaviral'] )

workflow test_spades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    SPADES ( input, [] )
}

workflow test_spades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]

    SPADES ( input, [] )
}
