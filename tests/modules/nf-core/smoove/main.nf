#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SMOOVE } from '../../../../modules/nf-core/smoove/main.nf'

workflow test_smoove {
    
    bam_tuple_ch = Channel.of([ [ id:'test' ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)
    				])
    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)
    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    SMOOVE ( bam_tuple_ch, fasta, fai )
}
