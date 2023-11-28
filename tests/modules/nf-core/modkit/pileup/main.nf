#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MODKIT_PILEUP } from '../../../../../modules/nf-core/modkit/pileup/main.nf'

workflow test_modkit_pileup {
    
    bam_input = [
        [ id:'test', single_end:false ],
        file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    bai_input = [
    file(params.test_data['sarscov2']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true) ]
    output_bed = [ file("test.bed") ]

    MODKIT_PILEUP ( bam_input, bai_input, output_bed )
}