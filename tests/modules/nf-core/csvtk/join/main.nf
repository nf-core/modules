#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CSVTK_JOIN } from '../../../../../modules/nf-core/csvtk/join/main.nf'

workflow test_csvtk_join {
    
    input = [ 
        [ id:'test' ], // meta map
        [ 
            file("https://github.com/nf-core/test-datasets/raw/bacass/bacass_hybrid.csv", checkIfExists: true),
            file("https://github.com/nf-core/test-datasets/raw/bacass/bacass_short.csv", checkIfExists: true) 
        ]
    ]

    CSVTK_JOIN ( input )
}
