#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HIFIASM } from '../../../modules/hifiasm/main.nf'

/*
 * Test with long reads only
 */
workflow test_hifiasm_hifi_only {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
    ]

    HIFIASM ( input, [], [], false )
}

/*
 * Test with parental reads for phasing
 */
workflow test_hifiasm_with_parental_reads {
    input = [
        [ id:'test' ], // meta map
        [ file(params.test_data['homo_sapiens']['pacbio']['hifi'], checkIfExists: true) ]
    ]
    paternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test_yak'], checkIfExists: true)
    maternal_kmer_dump = file(params.test_data['homo_sapiens']['illumina']['test2_yak'], checkIfExists: true)

    HIFIASM ( input, paternal_kmer_dump, maternal_kmer_dump, true )
}
