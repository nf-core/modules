#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EPANG } from "$moduleDir/modules/nf-core/epang/main.nf"

workflow test_epang {

    input = [
        [ id:'test', model:'LG' ], // meta map
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/query.alnfaa.gz', checkIfExists: true)
    ]

    EPANG (
        input,
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/reference.alnfaa.gz', checkIfExists: true),
        file('https://github.com/nf-core/test-datasets/raw/modules/data/delete_me/epang/reference.newick', checkIfExists: true),
        [],
        [],
        []
    )
}
