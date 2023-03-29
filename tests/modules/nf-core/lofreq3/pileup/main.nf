#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LOFREQ3_PILEUP } from '../../../../../modules/nf-core/lofreq3/pileup/main.nf'

workflow test_lofreq3_pileup {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    LOFREQ3_PILEUP ( input )
}
