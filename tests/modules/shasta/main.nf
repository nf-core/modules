#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHASTA } from '../../../modules/shasta/main.nf'

workflow test_shasta {
    
    input = [ 
        [ id:'test', model:'Nanopore-Oct2021' ], // meta map
        [ file(params.test_data['bacteroides_fragilis']['nanopore']['test_fastq_gz']) ],
    ]

    SHASTA ( input )
}
