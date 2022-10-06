#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KRAKENUNIQ_KRAKENUNIQ } from '../../../../../modules/nf-core/krakenuniq/krakenuniq/main.nf'

workflow test_krakenuniq_krakenuniq {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    KRAKENUNIQ_KRAKENUNIQ ( input )
}
