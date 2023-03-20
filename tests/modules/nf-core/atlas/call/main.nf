#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATLAS_CALL } from '../../../../../modules/nf-core/atlas/call/main.nf'

workflow test_atlas_call {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        [],
        []
    ]
    fasta         = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
    fai           = file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true)
    known_alleles = []
    method        = 'randomBase'

    ATLAS_CALL ( input, fasta, fai, known_alleles, method )
}
