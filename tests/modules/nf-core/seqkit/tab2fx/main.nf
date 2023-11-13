#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQKIT_TAB2FX                               } from '../../../../../modules/nf-core/seqkit/tab2fx/main.nf'
include { SEQKIT_TAB2FX as SEQKIT_TAB2FX_GZ           } from '../../../../../modules/nf-core/seqkit/tab2fx/main.nf'
include { SEQKIT_TAB2FX as SEQKIT_TAB2FX_FQ           } from '../../../../../modules/nf-core/seqkit/tab2fx/main.nf'
include { SEQKIT_TAB2FX as SEQKIT_TAB2FX_FQ_GZ        } from '../../../../../modules/nf-core/seqkit/tab2fx/main.nf'

workflow test_seqkit_tab2fx_fasta {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta_txt_zst'], checkIfExists: true) ]
            ]

    SEQKIT_TAB2FX ( input )

}

workflow test_seqkit_tab2fx_fastq {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_txt_zst'], checkIfExists: true) ]
            ]

    SEQKIT_TAB2FX_FQ ( input )

}

workflow test_seqkit_tab2fx_fasta_gz {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta_txt_zst'], checkIfExists: true) ]
            ]

    SEQKIT_TAB2FX_GZ ( input )

}

workflow test_seqkit_tab2fx_fastq_gz {

    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_txt_zst'], checkIfExists: true) ]
            ]

    SEQKIT_TAB2FX_FQ_GZ ( input )

}
