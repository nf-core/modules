#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_FX2TAB                            } from '../../../../../modules/nf-core/seqkit/fx2tab/main.nf'

workflow test_seqkit_fx2tab_fasta {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]

    SEQKIT_FX2TAB ( input )

}

workflow test_seqkit_fx2tab_fastq {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]

    SEQKIT_FX2TAB ( input )

}
