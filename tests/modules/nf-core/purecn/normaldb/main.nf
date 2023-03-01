#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_NORMALDB } from '../../../../../modules/nf-core/purecn/normaldb/main.nf'

workflow test_purecn_normaldb {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURECN_NORMALDB ( input )
}
