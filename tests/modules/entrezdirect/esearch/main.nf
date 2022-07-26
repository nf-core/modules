#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTREZDIRECT_ESEARCH } from '../../../../modules/entrezdirect/esearch/main.nf'

workflow test_entrezdirect_esearch {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ENTREZDIRECT_ESEARCH ( input )
}
