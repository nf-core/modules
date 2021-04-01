#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPADES } from '../../../software/spades/main.nf' addParams( spades_hmm: false ,options: ['args': '--rnaviral'] )

workflow test_spades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    coronaspades = false

    SPADES ( input, [], coronaspades )
}

workflow test_spades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    coronaspades = false

    SPADES ( input, [], coronaspades )
}

workflow test_coronospades_single_end {
    input = [ [ id:'test', single_end:true ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    coronaspades = true

    SPADES ( input, [], coronaspades )
}

workflow test_coronospades_paired_end {
    input = [ [ id:'test', single_end:false ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
            ]
    coronaspades = true

    SPADES ( input, [], coronaspades )
}

