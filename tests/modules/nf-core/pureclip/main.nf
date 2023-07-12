#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECLIP } from '../../../../modules/nf-core/pureclip/main.nf'

workflow test_pureclip {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURECLIP ( input )
}
