#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EIGENSTRATDATABASETOOLS_SNPCOVERAGE } from '../../../../modules/eigenstratdatabasetools/snpcoverage/main.nf'

workflow test_eigenstratdatabasetools_snpcoverage {

    input = [
        [ id:'test', single_end:false ], // meta map
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/pileupcaller.double.geno.txt", checkIfExists: true),
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/pileupcaller.double.snp.txt", checkIfExists: true),
        file("/Users/lamnidis/Software/github/TCLamnidis/modules/local_test_data/pileupcaller.double.ind.txt", checkIfExists: true)
    ]

    EIGENSTRATDATABASETOOLS_SNPCOVERAGE ( input )
}
