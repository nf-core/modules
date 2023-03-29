#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ3_VARCALLING } from '../../../../../modules/nf-core/lofreq3/varcalling/main.nf'

workflow test_lofreq3_varcalling {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    LOFREQ3_VARCALLING ( input )
}
