#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FQ_GENERATE } from '../../../../../modules/nf-core/fq/generate/main.nf'

workflow test_fq_generate {
    
    input = Channel.of(
        [ id:'test', single_end:false ]
    )

    FQ_GENERATE ( input )
}
