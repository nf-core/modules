#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARRIBA_DOWNLOAD } from '../../../../../modules/nf-core/arriba/download/main.nf'

workflow test_arriba_download {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ARRIBA_DOWNLOAD ( input )
}
