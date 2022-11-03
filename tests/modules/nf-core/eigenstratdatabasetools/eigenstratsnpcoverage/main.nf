#!/usr/bin/env nextflow

// Will need a way to create an eigenstrat dataset to test these tools on. Need pileupcaller module.
nextflow.enable.dsl = 2

include { SEQUENCETOOLS_PILEUPCALLER } from '../../../../../modules/nf-core/sequencetools/pileupcaller/main.nf'
include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE } from '../../../../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main.nf'

workflow test_eigenstratdatabasetools_eigenstratsnpcoverage {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_mpileup'], checkIfExists: true)
    ]

    snpfile = file(params.test_data['homo_sapiens']['genome']['genome_21_eigenstrat_snp'], checkIfExists: true)

    SEQUENCETOOLS_PILEUPCALLER ( input, snpfile, [] )
    EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE ( SEQUENCETOOLS_PILEUPCALLER.out.eigenstrat )
}
