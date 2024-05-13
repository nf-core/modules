#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EPANG_PLACE } from '../../../../../modules/nf-core/epang/place/main.nf'

workflow test_epang_place {
    
    input = [
        [ id:'test', model:'LG' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/query.alnfaa.gz', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/reference.alnfaa.gz', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/reference.newick', checkIfExists: true)
    ]

    EPANG_PLACE ( 
        input,
        [],
        []
    )
}
