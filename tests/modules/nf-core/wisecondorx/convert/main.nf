#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { WISECONDORX_CONVERT } from '../../../../../modules/nf-core/wisecondorx/convert/main.nf'

workflow test_wisecondorx_convert_bam {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true)
    ]

    WISECONDORX_CONVERT (
        input,
        [[],[]],
        [[],[]]
    )
}

workflow test_wisecondorx_convert_cram {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_cram_crai'], checkIfExists: true)
    ]

    fasta = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    ]

    fasta_fai = [
        [ id:'test' ],
        file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)
    ]

    WISECONDORX_CONVERT (
        input,
        fasta,
        fasta_fai
    )
}
