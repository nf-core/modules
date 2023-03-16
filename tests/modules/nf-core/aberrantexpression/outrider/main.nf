#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ABERRANTEXPRESSION_OUTRIDER } from '../../../../../modules/nf-core/aberrantexpression/outrider/main.nf'

workflow test_aberrantexpression_outrider {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    ABERRANTEXPRESSION_OUTRIDER ( input )
}
