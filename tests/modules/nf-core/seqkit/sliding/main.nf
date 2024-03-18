#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_SLIDING                            } from '../../../../../modules/nf-core/seqkit/sliding/main.nf'

workflow test_seqkit_sliding_fasta {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]

    SEQKIT_SLIDING ( input )

}

workflow test_seqkit_sliding_fastq {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    SEQKIT_SLIDING ( input )

}
