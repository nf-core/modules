#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHASTA } from '../../../../modules/nf-core/shasta/main.nf'

workflow test_shasta {
    
    input = [ 
        [ id:'test', model:'Nanopore-Dec2019' ], // meta map
        [ file(params.test_data['candidatus_portiera_aleyrodidarum']['nanopore']['test_fastq_gz']) ],
    ]

    SHASTA ( input )
}
