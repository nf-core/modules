#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PEPPER_MARGIN_DEEPVARIANT } from '../../../../../modules/nf-core/pepper/pepper_margin_deepvariant/main.nf'

workflow test_pepper_margin_deepvariant {

    bam_tuple_ch = Channel.of([ [ id:'test' ], // meta map
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam'], checkIfExists: true),
                               file(params.test_data['homo_sapiens']['illumina']['test_paired_end_sorted_bam_bai'], checkIfExists: true),
                               []
                                ])

    fasta = file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true)

    fai = file(params.test_data['homo_sapiens']['genome']['genome_fasta_fai'], checkIfExists: true)

    PEPPER_MARGIN_DEEPVARIANT ( bam_tuple_ch, fasta, fai )
}
