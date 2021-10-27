#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEQTK_MERGEPE } from '../../../../modules/seqtk/mergepe/main.nf' addParams( options: [:] )

workflow test_seqtk_mergepe {
    
    input = [ [ id:'test', single_end:false ], // meta map
              file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true) ]

    SEQTK_MERGEPE ( input )
}
