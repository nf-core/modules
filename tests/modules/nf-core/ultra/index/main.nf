#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRA_INDEX } from '../../../../../modules/nf-core/ultra/index/main.nf'
include { GNU_SORT    } from '../../../../../modules/nf-core/gnu/sort/main.nf'

workflow test_ultra_index {

    input = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
    ]
    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    GNU_SORT ( input )
    ULTRA_INDEX ( genome, GNU_SORT.out.sorted.map{ it[1] } )
}
