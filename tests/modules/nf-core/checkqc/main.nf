#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHECKQC } from '../../../../modules/nf-core/checkqc/main.nf'

workflow test_checkqc {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    CHECKQC ( input )
}
