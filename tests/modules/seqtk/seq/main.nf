#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ } from '../../../../modules/seqtk/seq/main.nf'

workflow test_seqtk_seq {
    
    input = [ 
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) 
    ]

    SEQTK_SEQ ( input )
}
