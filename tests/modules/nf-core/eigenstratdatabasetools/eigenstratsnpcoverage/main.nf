#!/usr/bin/env nextflow

// Will need a way to create an eigenstrat dataset to test these tools on. Need pileupcaller module.
nextflow.enable.dsl = 2

include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE } from '../../../../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main.nf'

workflow test_eigenstratdatabasetools_eigenstratsnpcoverage {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/Ancients+Saami.geno", checkIfExists: true),
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/Ancients+Saami.snp", checkIfExists: true),
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/Ancients+Saami.ind", checkIfExists: true)
    ]

    EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE ( input )
}
