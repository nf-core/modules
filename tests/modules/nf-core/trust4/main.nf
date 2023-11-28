#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRUST4 } from '../../../../modules/nf-core/trust4/main.nf'

workflow test_trust4 {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    TRUST4 ( input )
}
