#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NGMERGE } from '../../../../modules/nf-core/ngmerge/main.nf'

workflow test_ngmerge {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    NGMERGE ( input )
}
