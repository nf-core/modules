#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EKLIPSE } from '../../../../modules/nf-core/eklipse/main.nf'

workflow test_eklipse {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        '/usr/local/bin/data/test_illumina.bam'
    ]

    EKLIPSE ( input, [] )
}
