#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EMBOSS_SEQRET } from '../../../../modules/emboss/seqret/main.nf'

workflow test_emboss_seqret {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    EMBOSS_SEQRET ( input )
}
