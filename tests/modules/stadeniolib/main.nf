#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { STADENIOLIB } from '../../../modules/stadeniolib/main.nf'

workflow test_stadeniolib {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    STADENIOLIB ( input, file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true), file(params.test_data['sarscov2']['genome']['genome__fasta_fai'], checkIfExists: true))
}
