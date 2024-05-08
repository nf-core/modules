#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ULTRA_PIPELINE } from '../../../../../modules/nf-core/ultra/pipeline/main.nf'
include { GUNZIP         } from '../../../../../modules/nf-core/gunzip/main.nf'
include { GNU_SORT       } from '../../../../../modules/nf-core/gnu/sort/main.nf'

workflow test_ultra_pipeline {

    input_gunzip = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true)
    ]

    input_sort = [
        [ id:'test', single_end:false ],
        file(params.test_data['homo_sapiens']['genome']['genome_gtf']  , checkIfExists: true)
    ]

    genome = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    GUNZIP ( input_gunzip )
    GNU_SORT ( input_sort )
    ULTRA_PIPELINE (
        GUNZIP.out.gunzip,
        genome,
        GNU_SORT.out.sorted.map{ it[1] } )

}
