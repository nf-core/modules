#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_RENAME } from '../../../../modules/seqtk/rename/main.nf'

workflow test_seqtk_rename {
    input = [ [ id:'test' ],   // meta map
              [ file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true) ]
            ]
    SEQTK_RENAME ( input )
}

workflow test_seqtk_rename_fq {
    input = [ [ id:'test' ], // meta map
              [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
            ]
    SEQTK_RENAME ( input )
}
