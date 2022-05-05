#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SHASTA } from '../../../modules/shasta/main.nf'

workflow test_shasta {
    
    input = [ 
        [ id:'test', model:'Nanopore-Oct2021' ], // meta map
        [ file('https://github.com/nf-core/test-datasets/raw/bacass/nanopore/subset15000.fq.gz') ],
    ]

    SHASTA ( input )
}
