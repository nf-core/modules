#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MAPAD_INDEX } from '../../../../../modules/nf-core/mapad/index/main.nf'
include { MAPAD_MAP } from '../../../../../modules/nf-core/mapad/map/main.nf'

//
// Test with single-end data
//
workflow test_mapad_map_single_end {
    input = [
        [ id:'test', single_end:true ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_single_end_bam'], checkIfExists: true)
        ]
    ]
    fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)

    mismatch_parameter = 0.03
    double_stranded_library = false
    five_prime_overhang = 0.5
    three_prime_overhang = 0.5
    deam_rate_double_stranded = 0.02
    deam_rate_single_stranded = 1.0
    indel_rate = 0.001

    MAPAD_INDEX ( [ [:], fasta ] )
    MAPAD_MAP ( input, MAPAD_INDEX.out.index, mismatch_parameter, double_stranded_library, five_prime_overhang, three_prime_overhang, deam_rate_double_stranded, deam_rate_single_stranded, indel_rate )
}