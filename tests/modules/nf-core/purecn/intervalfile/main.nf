#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PURECN_INTERVALFILE } from '../../../../../modules/nf-core/purecn/intervalfile/main.nf'

workflow test_purecn_intervalfile {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    PURECN_INTERVALFILE ( input )
}
