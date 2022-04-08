#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SNAPALIGNER_PAIRED } from '../../../../modules/snapaligner/paired/main.nf'

workflow test_snapaligner_paired {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SNAPALIGNER_PAIRED ( input )
}
