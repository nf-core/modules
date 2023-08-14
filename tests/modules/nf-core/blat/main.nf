#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_SEQ } from '../../../../modules/nf-core/seqtk/seq/main.nf'
include { BLAT      } from '../../../../modules/nf-core/blat/main.nf'

workflow test_blat {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]
    sequences = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    
    SEQTK_SEQ ( input )
    BLAT ( SEQTK_SEQ.out.fastx, sequences )
}
