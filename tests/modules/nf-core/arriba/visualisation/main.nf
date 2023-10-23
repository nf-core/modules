#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ARRIBA_VISUALISATION } from '../../../../../modules/nf-core/arriba/visualisation/main.nf'

workflow test_arriba_visualisation {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ARRIBA_VISUALISATION ( input )
}
