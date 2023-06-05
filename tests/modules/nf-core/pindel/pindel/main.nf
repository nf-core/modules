#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PINDEL_PINDEL } from '../../../../../modules/nf-core/pindel/pindel/main.nf'

workflow test_pindel_pindel {
    
    input = [
        [ id:'test', single_end:false, insert_size:500 ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]

    PINDEL_PINDEL ( input, file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true), file(params.test_data['sarscov2']['genome']['genome_fasta_fai'], checkIfExists: true), file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true) )
}
