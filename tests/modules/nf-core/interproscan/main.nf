#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INTERPROSCAN } from '../../../../modules/nf-core/interproscan/main.nf'

workflow test_interproscan {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    INTERPROSCAN ( input )
}
