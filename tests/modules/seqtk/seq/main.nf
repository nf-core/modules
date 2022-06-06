#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ } from '../../../../modules/seqtk/seq/main.nf'

workflow test_seqtk_seq {
    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]
    SEQTK_SEQ ( input )
}

workflow test_seqtk_seq_fq {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    SEQTK_SEQ ( input )
}