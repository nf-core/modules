#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SENTIEON_TNSCOPE              } from '../../../../../modules/nf-core/sentieon/tnscope/main.nf'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../../../modules/nf-core/picard/addorreplacereadgroups/main.nf'

workflow test_tnscope {

    fasta = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)]).collect()
    fai   = Channel.of([ [ id:'test' ], file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)]).collect()
    bai   = Channel.of([file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true)])

    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true)
    ]
    PICARD_ADDORREPLACEREADGROUPS ( input )

    ch_tnscope_in = PICARD_ADDORREPLACEREADGROUPS.out.bam.combine(bai)

    SENTIEON_TNSCOPE ( ch_tnscope_in, fasta, fai, [[:],[],[]], [[:],[],[]], [[:],[],[]], [[:],[]] )
}
