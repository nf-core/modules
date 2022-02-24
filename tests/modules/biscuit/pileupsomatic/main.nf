#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BISCUIT_INDEX         } from '../../../../modules/biscuit/index/main.nf'
include { BISCUIT_PILEUPSOMATIC } from '../../../../modules/biscuit/pileupsomatic/main.nf'

workflow test_biscuit_pileupsomatic {

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
        file(params.test_data['homo_sapiens']['illumina']['test2_paired_end_sorted_bam_bai'], checkIfExists: true)
    ]
    fasta   = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    BISCUIT_INDEX( fasta )
    BISCUIT_PILEUPSOMATIC ( input, BISCUIT_INDEX.out.index )
}
