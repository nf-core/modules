#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMNCOPYNUMBERCALLER } from '../../../../modules/nf-core/smncopynumbercaller/main.nf'

workflow test_smncopynumbercaller {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    SMNCOPYNUMBERCALLER ( input )
}
