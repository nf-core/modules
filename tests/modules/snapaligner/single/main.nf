#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNAPALIGNER_SINGLE } from '../../../../modules/snapaligner/single/main.nf'

workflow test_snapaligner_single {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SNAPALIGNER_SINGLE ( input )
}
